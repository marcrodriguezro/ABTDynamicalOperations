#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <chrono>
#include "abt_compression.h"

using namespace std;

bool Input(const char *filename, vector<pair<int, int>> *edges);

int main(/*int argc, char **argv*/)
{
  /*if (argc != 3)
  {
    cerr << "Incorrect number of arguments check usage for info" << endl;
    exit(EXIT_FAILURE);
  }*/
    
  vector<int> str;
  vector<pair<int, int>> edgesReconstructed;
  vector<pair<int, int>> edges;
  if (!Input("../data/graph/modified/tvshow_edges3.txt", &edges))
  {
    cerr << "error: Load failed" << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Edges      : " << edges.size() << endl;

 
  AbtCompression bl;
  string result="test.txt";
  int val1_1, val1_2;
  auto start = std::chrono::high_resolution_clock::now();
  bl.Compress(edges, result);

  int start_index, n, height;
  std::sort(edges.begin(), edges.end());
  auto max_pair = *std::max_element(std::begin(edges), std::end(edges), [](const auto& p1, const auto& p2) {
      return std::max(p1.first, p1.second) < std::max(p2.first, p2.second);
      });
  int max = std::max(max_pair.first, max_pair.second) + 1;
  height = (int)ceil(log2(max)) + 1;
  n = std::pow(2, height) - 1;


  streampos begin, end;
  ifstream myfile(result.c_str(), ios::binary);
  begin = myfile.tellg();
  //myfile.seekg(0, ios::end);
  end = myfile.tellg();
  myfile.close();

  //bl.DecodeBit(result, str);

 
  long size = (end - begin) * 8;
  cout << "Bit length : " << size << endl;
  cout << "Bits/edge  : " << ((double)size / edges.size()) << endl;
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  cout << "Elapsed time to compress: " << elapsed.count() << " s\n";
  //-------------------------------------------------------------------
  cout << "Enter a query to search for. (First enter a number, followed by an space and then enter the second number " << endl;
  cin >> val1_1 >> val1_2;
  //bl.Dynamic_decoder(height, str, &edgesReconstructed);
  bl.Partial_decompression(result, val1_1, val1_2, n);
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
