#include "abt_compression.h"
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include<list>
#include <sstream>
#include <string.h>

using namespace std;

class bitChar
{
public:
  unsigned char *c;
  int shift_count;
  string BITS;

  bitChar()
  {
    shift_count = 0;
    c = (unsigned char *)calloc(1, sizeof(char));
  }

  string readByBits(ifstream &inf)
  {
    string s = "";
    char buffer[1];
    while (inf.read(buffer, 1))
    {
      s += getBits(*buffer);
    }
    return s;
  }

  void setBITS(string X)
  {
    BITS = X;
  }

  int insertBits(ofstream &outf)
  {
    int total = 0;

    while (BITS.length())
    {
        if (BITS[0] == '1')
            *c |= 1;
        *c <<= 1;
        ++shift_count;
        ++total;
        BITS.erase(0, 1);

        if (shift_count == 7)
        {
            if (BITS.size() > 0)
            {
                if (BITS[0] == '1')
                    *c |= 1;
                ++total;
                BITS.erase(0, 1);
            }

            writeBits(outf);
            shift_count = 0;
            free(c);
            c = (unsigned char*)calloc(1, sizeof(char));
        }
    }

    if (shift_count > 0)
    {
      *c <<= (7 - shift_count);
      writeBits(outf);
      free(c);
      c = (unsigned char *)calloc(1, sizeof(char));
    }
    outf.close();
    return total;
  }

  string getBits(unsigned char X)
  {
    stringstream itoa;
    for (unsigned s = 7; s > 0; s--)
    {
      itoa << ((X >> s) & 1);
    }

    itoa << (X & 1);
    return itoa.str();
  }

  void writeBits(ofstream &outf)
  {
    outf << *c;
  }

  ~bitChar()
  {
    if (c)
      free(c);
  }
};

bool BitString::Input(const char *filename)
{
  std::ifstream ifs(filename, std::ios::binary);
  if (!ifs)
    return false;
  data_.clear();
  uint64_t d;
  while (ifs.read((char *)&d, sizeof(uint64_t)))
  {
    data_.push_back(d);
  }
  length_ = (uint64_t)data_.size() * 64;
  return !ifs.bad();
}

bool BitString::Output(const char *filename)
{
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs)
    return false;
  for (uint64_t i = 0; i < length_; i += 64)
  {
    ofs.write((char *)&data_[i >> 6], sizeof(uint64_t));
  }
  return ofs.good();
}

void DeltaCode::EncodeInt(int val, BitString *out)
{
  out->Init(0);
  int length_val = CountBitLength(val + 1) - 1;
  int length_length_val = CountBitLength(length_val + 1) - 1;
  for (int i = 0; i < length_length_val; ++i)
    out->AppendBit(0);
  out->AppendBit(1);
  for (int i = length_length_val - 1; i >= 0; --i)
  {
    out->AppendBit((length_val + 1) >> i & 1);
  }
  for (int i = length_val - 1; i >= 0; --i)
  {
    out->AppendBit((val + 1) >> i & 1);
  }
}

int DeltaCode::DecodeNextInt(const BitString &in, uint64_t *i)
{
  int length_length_val = 0;
  while (!in.GetBit(*i))
  {
    ++length_length_val;
    ++(*i);
  }
  ++(*i);
  int length_val = 1;
  for (int j = 0; j < length_length_val; ++j)
  {
    length_val = (length_val << 1) | in.GetBit(*i);
    ++(*i);
  }
  --length_val;
  int val = 1;
  for (int j = 0; j < length_val; ++j)
  {
    val = (val << 1) | in.GetBit(*i);
    ++(*i);
  }
  return val - 1;
}

void AbtCompression ::TransformToAdj(const std::vector<std::pair<int, int>> &edges,
                                     bool directed, std::vector<std::vector<int>> *adj)
{
  int num_v = 0;
  //std::cout << "\nEdges.size() : " << edges.size();
  for (uint64_t i = 0; i < edges.size(); ++i)
  {
    num_v = std::max(num_v, std::max(edges[i].first, edges[i].second) + 1);
  }
  //std::cout << "\nNum_v : " << num_v;
  adj->resize(num_v);
  for (uint64_t i = 0; i < edges.size(); ++i)
  {
    adj->at(edges[i].first).at(edges[i].second) = 1;
    //std::cout << std::endl
    //<< "adj->at(edges[i].first) : " << &adj->at(edges[i].first);
    //std::cout << std::endl
    //<< "edges[i].second : " << edges[i].second;
  }
  if (!directed)
  {
    std::cout << std::endl
              << "enters";
    for (uint64_t i = 0; i < edges.size(); ++i)
    {
      adj->at(edges[i].second).push_back(edges[i].first);
    }
  }
  /* for (int i = 0; i < num_v; ++i)
  {
    std::sort(adj->at(i).begin(), adj->at(i).end());
  } */
}

void AbtCompression ::TransformToEdge(const std::vector<std::vector<int>> &adj,
                                      std::vector<std::pair<int, int>> *edges)
{
  for (int i = 0; i < adj.size(); ++i)
  {
    for (int j = 0; j < adj[i].size(); ++j)
    {
      edges->push_back(std::make_pair(i, adj[i][j]));
    }
  }
}

inline int parent(int index)
{
  return (int)((index - 1) / 2);
}

inline int sibling(int index)
{
  if ((index % 2) == 1)
  {
    return index + 1;
  }
  else
  {
    return index - 1;
  }
}
void AbtCompression ::Compress(std::vector<std::pair<int, int>> edges, string result)
{
  int current_index, start_index, n, height;
  bool dcn_reached;

  std::vector<int> order;
  std::sort(edges.begin(), edges.end());
  auto max_pair = *std::max_element(std::begin(edges), std::end(edges), [](const auto &p1, const auto &p2) {
    return std::max(p1.first, p1.second) < std::max(p2.first, p2.second);
  });
  int max = std::max(max_pair.first, max_pair.second) + 1;
  height = (int)ceil(log2(max)) + 1;
  n = std::pow(2, height) - 1;
  start_index = (int)n / 2;

  std::vector<std::vector<int>>
  input_array(max, std::vector<int>(max, 0));
  TransformToAdj(edges, true, &input_array);

  for (int i = 0; i < max; ++i)
  {
    std::vector<int> temp_array(n, -1);
    for (int j = 0; j < max; ++j)
    {
      if (input_array[i][j] == 1)
      {
        current_index = start_index + j;
        dcn_reached = false;
        while (dcn_reached == false)
        {
          temp_array[current_index] = 1;
          if (temp_array[parent(current_index)] != 1)
          {
            temp_array[parent(current_index)] = 1;
            temp_array[sibling(current_index)] = 0;
            current_index = parent(current_index);
          }
          else
          {
            dcn_reached = true;
          }
        }
      }
    }

    ofstream outf(result.c_str(), std::ios_base::app);
    std::string str;
    for (int k = 0; k < n; k++)
    {
      if (temp_array[k] != -1)
      {
        str.push_back('0' + temp_array[k]);
      }
    }
    bitChar bchar;
    bchar.setBITS(str);
    bchar.insertBits(outf);
  }
}


void AbtCompression::Partial_decompression(std::string result, int u, int v, int n) {
    ifstream inf(result.c_str(), std::ios_base::app);
    bitChar bchar;
    std::string decoded;
    decoded = bchar.readByBits(inf);
    std::string aux;
   
    int cur, begCol, endCol, nextIndex, nodesAtDepth, nodesAtNextDepth;
    cur = 0;
    begCol = 0;
    endCol = n;
    nextIndex = 0;
    nodesAtDepth = 1;
    int midCol;
    nodesAtNextDepth = 0;
    int curNode;
    int pos;
    for (int i = 0; i < decoded.size();i= i +nodesAtDepth) {
        for (int j = i; j < nodesAtDepth; j++) {
            aux = decoded.substr(j, 1);
            curNode = stoi(aux, nullptr, 2);
            if (j == nextIndex) {
                midCol = (begCol + endCol) / 2;
                if (curNode==1) {
                    if (v < midCol) {
                        pos = 0;
                    }
                    else {
                        pos = 1;
                    }
                    nextIndex = i + nodesAtDepth + nodesAtNextDepth + pos;
                    if (pos == 1) {
                        begCol = midCol;
                    }
                    else {
                        endCol = midCol;
                    }
                    nodesAtNextDepth = nodesAtNextDepth + 2;
                }
                else if(curNode==0){
                    //result.at(j,1); // metodo para a�adir al bitstring el valor 1;
                    cout << "A�ado un nodo" << endl;
                    return;
                }
            }
            else if (curNode == 1) {
                nodesAtNextDepth = nodesAtNextDepth + 2;
            }
        }
        nodesAtDepth = nodesAtNextDepth;
        nodesAtNextDepth = 0;
    }

}



/*
void AbtCompression ::DecodeBit(string result, std::vector<int> &str) {
    ifstream inf(result.c_str(), std::ios_base::app);
    bitChar bchar;
    std::string decoded;
    decoded = bchar.readByBits(inf);
    int n = decoded.length();
    int *res = new int[n];
    int prueba;
    std::string aux;
    for (int i = 0; i < n; i++) { //int i = n; i >= 0; i--

        aux = decoded.substr(i, 1);
        if (aux != "") {
            str.push_back(stoi(aux, nullptr,2));
        }
    }
}*/
void AbtCompression ::Develop(const BitString &code, std::vector<std::pair<int, int>> *edges)
{
  DeltaCode delta;
  uint64_t cur = 0;
  int num_v = delta.DecodeNextInt(code, &cur);
  std::vector<int> order(num_v);
  for (int i = 0; i < num_v; ++i)
  {
    order[delta.DecodeNextInt(code, &cur)] = i;
  }
  std::vector<std::vector<int>> adj(num_v);
  for (int i = 0; i < num_v; ++i)
  {
    int j = i - delta.DecodeNextInt(code, &cur);
    if (code.GetBit(cur++))
      adj[i].push_back(i);
    // copying
    if (i != j)
    {
      for (int k = 0; k < adj[j].size(); ++k)
      {
        if (Exist(adj, adj[j][k], j) && j >= adj[j][k])
          continue;
        if (code.GetBit(cur++))
        {
          adj[i].push_back(adj[j][k]);
        }
      }
    }
    // residual
    int num_residual = delta.DecodeNextInt(code, &cur);
    if (num_residual != 0)
    {
      int sign = (code.GetBit(cur++) ? 1 : -1);
      int now = i + delta.DecodeNextInt(code, &cur) * sign;
      adj[i].push_back(now);
      for (int j = 1; j < num_residual; ++j)
      {
        now += delta.DecodeNextInt(code, &cur);
        adj[i].push_back(now);
      }
    }
    std::sort(adj[i].begin(), adj[i].end());
    // reciprocal
    for (int k = 0; k < adj[i].size(); ++k)
    {
      if (i >= adj[i][k])
        continue;
      if (code.GetBit(cur++))
      {
        adj[adj[i][k]].push_back(i);
      }
    }
  }
  TransformToEdge(adj, edges);
  for (uint64_t i = 0; i < edges->size(); ++i)
  {
    edges->at(i).first = order[edges->at(i).first];
    edges->at(i).second = order[edges->at(i).second];
  }
  std::sort(edges->begin(), edges->end());
}
/*
inline int recursive_Search(int i, std::vector<int>& str, std::vector<std::pair<int, int>>* edgesReconstructed, int height) {
    bool visited_left = false;
    bool visited_right = false;
    std::vector< int > visited_nodes;
    std::vector<int>::iterator it;
    int n = std::pow(2, height) - 1;
    int altura = 0;
    int changed_altura = 0;
    while (i != 14) {
        if (str[i] == 1 && i < (n-1) &&((i*2)+1) < n && visited_left == false) {

            it = std::find(visited_nodes.begin(), visited_nodes.end(), i);
            if (visited_nodes.empty()) {
                visited_nodes.push_back(0);
            }
            else if (*it != visited_nodes.back()) {
                visited_nodes.push_back(i);
            }
            i = (i * 2) + 1; //left sibling
            visited_nodes.push_back(i);
            altura++;
            //end para mirar que no sea null, back para el ultimo
            visited_right = false;
           /* if (str[i] == NULL) {
                edgesReconstructed->push_back(make_pair(0, i)); // esto hay que cambiarlo
            }*/
   /*
        }else if (str[i] == 1 && ((i * 2) + 2)< n && visited_right == false) { // ahora mismo falla pq llega al 8, vuelve arriba y como
            it = std::find(visited_nodes.begin(), visited_nodes.end(), i);
            if (it == visited_nodes.end()) { // check if the iterator points to null
                visited_nodes.push_back(i);
            }
            i = (i * 2) + 2; // right sibling
            cout << i << endl;
            altura++;
            if (str[i] == 0) {
                visited_left = true;
                visited_right = true;
            }
            visited_left = false;
        }else { //no tiene hijos ese nodo y no hemos visitado a su hermano

            it = std::find(visited_nodes.begin(), visited_nodes.end(), i);
            if ( it == visited_nodes.end() ) {
                visited_left = true;
                visited_right = true;
                visited_nodes.push_back(i);
            }
            i = (i - 1) / 2; // padre del actual
            altura--;

        }
    }
    return i;
}
*/
/*
void AbtCompression::Dynamic_decoder(int height, std::vector<int>& str, std::vector<std::pair<int, int>> *edgesReconstructed) { //MINE
    
    int i = 0;
    recursive_Search(i, str, edgesReconstructed, height);
}
*/



void AbtCompression ::BFSOrder(const std::vector<std::vector<int>> &adj, std::vector<int> *order)
{
  int num_v = adj.size();
  std::queue<int> q;
  int k = 0;
  for (int i = 0; i < num_v; ++i)
  {
    if (order->at(i) != -1)
      continue;
    order->at(i) = k++;
    q.push(i);
    while (!q.empty())
    {
      int v = q.front();
      q.pop();
      for (int j = 0; j < adj[v].size(); ++j)
      {
        int u = adj[v][j];
        if (order->at(u) != -1)
          continue;
        order->at(u) = k++;
        q.push(u);
      }
    }
  }
}

void AbtCompression ::Order(const std::vector<std::pair<int, int>> &edges,
                            std::vector<int> *order)
{
  std::vector<std::vector<int>> adj;
  TransformToAdj(edges, false, &adj);
  int num_v = adj.size();
  order->resize(num_v);
  std::fill(order->begin(), order->end(), -1);
  switch (kORDERING)
  {
  case BFS:
  {
    BFSOrder(adj, order);
    break;
  }
  default:
  {
    for (int i = 0; i < num_v; ++i)
      order->at(i) = i;
  }
  }
}

void AbtCompression ::CompressVertexes(const std::vector<std::vector<int>> &adj,
                                       BitString *result)
{
  int num_v = adj.size();
  DeltaCode delta;
  BitString tmp;
  for (int i = 0; i < num_v; ++i)
  {
    BitString best;
    best.Init(0);
    for (int j = i; j > i - kWINDOW_WIDTH && j >= 0; --j)
    {
      BitString now;
      now.Init(0);
      delta.EncodeInt(i - j, &tmp);
      now.AppendBitString(tmp);
      if (Exist(adj, i, i))
      { // self-loop
        now.AppendBit(1);
      }
      else
      {
        now.AppendBit(0);
      }
      std::vector<int> residual;
      // copying
      int cur = 0;
      if (i != j)
      {
        for (int k = 0; k < adj[j].size(); ++k)
        {
          if (Exist(adj, adj[j][k], j) && j >= adj[j][k])
            continue;
          ProceedCopying(adj, adj[j][k], i, &cur, &residual);
          if (cur < adj[i].size() && adj[i][cur] == adj[j][k])
          {
            now.AppendBit(1);
          }
          else
          {
            now.AppendBit(0);
          }
        }
      }
      ProceedCopying(adj, num_v, i, &cur, &residual);
      // residual
      delta.EncodeInt(residual.size(), &tmp);
      now.AppendBitString(tmp);
      if (residual.size() != 0)
      {
        if (residual[0] > i)
        {
          now.AppendBit(1);
        }
        else
        {
          now.AppendBit(0);
        }
        delta.EncodeInt(std::abs(residual[0] - i), &tmp);
        now.AppendBitString(tmp);
        for (int k = 1; k < residual.size(); ++k)
        {
          delta.EncodeInt(residual[k] - residual[k - 1], &tmp);
          now.AppendBitString(tmp);
        }
      }
      // reciprocal
      for (int k = 0; k < adj[i].size(); ++k)
      {
        if (i >= adj[i][k])
          continue;
        if (Exist(adj, adj[i][k], i))
        {
          now.AppendBit(1);
        }
        else
        {
          now.AppendBit(0);
        }
      }
      // update best
      if (i == j || best.get_length() > now.get_length())
        best = now;
    }
    result->AppendBitString(best);
  }
}
