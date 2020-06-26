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
    int entered_value1 = 0;
    int entered_value2 = 0;

    cout << "If you want to check the existance of a node press 1, otherwise if you want to add a node press 2, if you want to check the existance of a node with a non-dynamical algorithm press 3, finally if you want to add a node in a static way press 4" << endl;
    cin >> button_pressed;
    cout << "If you want to stop the execution enter -1" << endl;
    if (button_pressed == 1) {
        while (entered_value1 != -1 || entered_value2 != -1) {
            cout << "Enter the leaf nodes that you want to check if exists any conection between them on the binary tree" << endl;
            cin >> entered_value1;
            cin >> entered_value2;
            if (entered_value1 > n / 2 || entered_value2 > n / 2) {
                cout << "Some of the nodes that you were searching DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value1 < 0 || entered_value2 < 0) {
                cout << "Some of the nodes that you were searching DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                auto start_dynamic = std::chrono::high_resolution_clock::now();
                check_existance_v = bl.Node_existance_checking(result, entered_value1, entered_value2, height);
                if (check_existance_v) {
                    cout << "The nodes and the conection between those nodes that you were searching EXIST" << endl;
                    auto finish_dynamically = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_dynamically - start_dynamic;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
                else {
                    cout << "Seems like the conection between the nodes that you entered DON'T EXIST" << endl;
                    auto finish_dynamically = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_dynamically - start_dynamic;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
            }
        }
    }
    else if(button_pressed == 2){
        while (entered_value1 != -1 || entered_value2 != -1) {
            cout << "Enter the conection between the nodes that you want to add" << endl;
            cin >> entered_value1;
            cin >> entered_value2;
            if (entered_value1 > n / 2 || entered_value2 > n / 2) {
                cout << "Some of the nodes that you are refering DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value1 < 0 || entered_value2 < 0) {
                cout << "Some of the nodes that you are refering DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                auto start_dynamically = std::chrono::high_resolution_clock::now();
                check_existance_v = bl.node_addition(result, entered_value1, entered_value2, height, max);

                if (check_existance_v) {
                    cout << "The the conection between those nodes has been added correclty." << endl;
                    auto finish_dynamically = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_dynamically - start_dynamically;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
                else {
                    cout << "There's some trouble adding the nodes and the conection between those nodes you selected" << endl;
                    auto finish_dynamically = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_dynamically - start_dynamically;
                    cout << "Elapsed time to add an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
            }
        }

    }
    else if (button_pressed == 3) {
        while (entered_value1 != -1 || entered_value2 != -1) {
            cout << "Enter the leaf nodes that you want to check if exists any conection between them on the binary tree" << endl;
            cin >> entered_value1;
            cin >> entered_value2;
            if (entered_value1 > n / 2 || entered_value2 > n / 2) {
                cout << "Some of the nodes that you were searching DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value1 < 0 || entered_value2 < 0) {
                cout << "Some of the nodes that you were searching DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                auto start_static = std::chrono::high_resolution_clock::now();
                check_existance_v = bl.Node_existance_checking_NCompress(result, entered_value1, entered_value2, height);
                if (check_existance_v) {
                    cout << "The nodes and the conection between those nodes that you were searching EXIST" << endl;
                    auto finish_static = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_static - start_static;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
                else {
                    cout << "Seems like the conection between the nodes that you entered DON'T EXIST" << endl;
                    auto finish_static = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_static - start_static;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
            }
        }
    }
    else if (button_pressed == 4) {
        while (entered_value1 != -1 || entered_value2 != -1) {
            cout << "Enter the conection between the nodes that you want to add in a static way" << endl;
            cin >> entered_value1;
            cin >> entered_value2;
            if (entered_value1 > n / 2 || entered_value2 > n / 2) {
                cout << "Some of the nodes that you are refering DO NOT EXIST because the tree is not that bigger " << endl;
            }
            else if (entered_value1 < 0 || entered_value2 < 0) {
                cout << "Some of the nodes that you are refering DO NOT EXIST because you entered a number lower than 0" << endl;
            }
            else {
                auto start_static = std::chrono::high_resolution_clock::now();
                check_existance_v = bl.node_addition(result, entered_value1, entered_value2, height, max);

                if (check_existance_v) {
                    cout << "The the conection between those nodes has been added correclty." << endl;
                    auto finish_static = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish_static - start_static;
                    cout << "Elapsed time to check the existance of an edge and a node dynamically: " << elapsed.count() << " s\n";
                }
                else {
                    cout << "There's some trouble adding the nodes and the conection between those nodes you selected" << endl;
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
