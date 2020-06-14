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
    std::ofstream ofs;
    ofs.open(result, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    //free what test.txt contains. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    int height;
    std::sort(edges.begin(), edges.end());
    auto max_pair = *std::max_element(std::begin(edges), std::end(edges), [](const auto& p1, const auto& p2) {
        return std::max(p1.first, p1.second) < std::max(p2.first, p2.second);
        });
    int max = std::max(max_pair.first, max_pair.second) + 1;
    height = (int)ceil(log2(max)) + 1;
    int n = std::pow(2, height) - 1; //num of nodes that the tree contains
    bool check_existance_v = false;
    int button_pressed;
    int entered_value = 0;

    cout << "If you want to check the existance of a node press 1, otherwise if you want to add a node press another number diferent than 1, for example: 2" << endl;
    cin >> button_pressed;
    cout << "If you want to stop the execution enter -1" << endl;
    if (button_pressed == 1) {
        while (entered_value != -1) {
            cout << "Enter the leaf node you want to check if exists on the binary tree" << endl;
            cin >> entered_value;
            cout << n <<endl;
            if (entered_value > n / 2) {
                cout << "The node you were searching DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value < 0) {
                cout << "The node you were searching DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                check_existance_v = bl.Node_existance_checking(result, entered_value, height);
                if (check_existance_v) {
                    cout << "The leaf node you were searching EXIST" << endl;
                }
                else {
                    cout << "The leaf node you were searching DO NOT EXIST" << endl;
                }
            }
        }
    }
    else {
        while (entered_value != -1) {
            cout << "Enter the leaf node you want to add on the binary tree" << endl;
            cin >> entered_value;
            if (entered_value > n / 2) {
                cout << "The node you are refering DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value < 0) {
                cout << "The node you are refering DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                check_existance_v = bl.node_addition(result, entered_value, height, max);
                if (check_existance_v) {
                    cout << "The node has been added correclty." << endl;
                }
                else {
                    cout << "There's some trouble adding the node you selected" << endl;
                }
            }
        }

    }
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
