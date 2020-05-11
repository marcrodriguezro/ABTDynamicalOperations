#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono>
#include "abt_compression.h"

using namespace std;

bool Input(const char *filename, vector<pair<int, int>> *edges);

int main(/*int argc, char **argv*/)
{
 /* if (argc != 3)
  {
    cerr << "Incorrect number of arguments check usage for info" << endl;
    exit(EXIT_FAILURE);
  }
  */
    std::vector<int> tmp;
  vector<pair<int, int>> edges;
  if (!Input("../data/graph/modified/tvshow_edges3.txt", &edges))
  {
    cerr << "error: Load failed" << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Edges      : " << edges.size() << endl;

  AbtCompression bl;
  string result = "test.txt";
  auto start = std::chrono::high_resolution_clock::now();
  bl.Compress(edges, result, tmp);
  streampos begin, end;
  ifstream myfile(result.c_str(), ios::binary);
  begin = myfile.tellg();
  myfile.seekg(0, ios::end);
  end = myfile.tellg();
  myfile.close();
  long size = (end - begin) * 8;
  cout << "Bit length : " << size << endl;
  cout << "Bits/edge  : " << ((double)size / edges.size()) << endl;
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  cout << "Elapsed time to compress: " << elapsed.count() << " s\n";

  int start_index, n, height;
  std::sort(edges.begin(), edges.end());
  auto max_pair = *std::max_element(std::begin(edges), std::end(edges), [](const auto& p1, const auto& p2) {
      return std::max(p1.first, p1.second) < std::max(p2.first, p2.second);
      });
  int max = std::max(max_pair.first, max_pair.second) + 1;
  height = (int)ceil(log2(max)) + 1;
  n = std::pow(2, height) - 1;

  bool check_existance_v = false;

  check_existance_v = bl.Edge_existance_checking(result, 7, height);
  //bl.Partial_decompression(result, 0, 7, n, tmp);
}

bool Input(const char *filename, vector<pair<int, int>> *edges)
{
  ifstream ifs(filename);
  if (!ifs)
    return false;
  for (int v, w; ifs >> v >> w;)
    edges->push_back(make_pair(v, w));
  return !ifs.bad();
}
