#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>


using std::vector;
using std::string;
using std::cout;
using std::endl;



uint32_t alphabetToNumerical(char alphabet){
    switch (alphabet){
        case '0': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        default: abort ();
    }
}

void readFasta(const string &filename, string &T){
    
    std::ifstream infile;
    infile.open(filename, std::ios::in);

    string line;
    while (getline(infile, line))
        if (line[0] != '>')
            for (char c : line)
                T.push_back(c);
    infile.close();
}

void radixSort(vector<uint32_t> &seq, vector<uint32_t> &ref, int n, int digit, int alphabet_size){
    
    vector<uint32_t> tmp(n, 0);
    // Sort the sequence according to the digit
    for (int j = digit-1; j >= 0; j--){
        // Count the number of each alphabet
        vector<uint32_t> count(alphabet_size+1, 0);
        for (int i = 0; i < n; i++)
            count[ref[seq[i] + j]]++;

        // Calculate the cummalative count of each element
        for (int i = 1; i < alphabet_size+1; i++)
            count[i] += count[i - 1];

        // Find the index of each element of the original array in count array, and place the elements in output array
        for (int i = n - 1; i >= 0; i--){
            tmp[count[ref[seq[i] + j]] - 1] = seq[i];
            count[ref[seq[i] + j]]--;
        }

        // Copy the sorted elements into original array
        for (int i = 0; i < n; i++)
            seq[i] = tmp[i];
    }
}

bool leq_pair(int a1, int a2, int b1, int b2) {
  return (a1 < b1 || (a1 == b1 && a2 <= b2));
}


bool leq_triple(int a1, int a2, int a3, int b1, int b2, int b3) {
  return (a1 < b1 || (a1 == b1 && leq_pair(a2, a3, b2, b3)));
}

void skewAlgorithm(vector<uint32_t> &s, vector<uint32_t> &SA, int n, int alphabet_size){
    
    // Calculate the length of mod suffixes
    int n0 = (n + 2)/3;
    int n1 = (n + 1)/3;
    int n2 = n/3;
    int n02 = n0 + n2;

    // Add 000 at the end
    s.push_back(0);
    s.push_back(0);
    s.push_back(0);

    /// Get SA12
    // Create triple indexes for sequence
    vector<uint32_t> s12(n02 + 3, 0);
    for (int i = 0, j = 0; i < n + (n0-n1); i++)
        if (i % 3 != 0)
            s12[j++] = i;
    
    // Sort the triplets based on letters
    vector<uint32_t> SA12(s12.begin(), s12.end());
    radixSort(SA12, s, n02, 3, alphabet_size);

    /// Get s12
    int name = 0;
    int c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++) {
        // Name the triplets based on the rank of triplet
        if (s[SA12[i]] != c0 || s[SA12[i] + 1] != c1 || s[SA12[i] + 2] != c2) {
            name++;
            c0 = s[SA12[i]];
            c1 = s[SA12[i] + 1];
            c2 = s[SA12[i] + 2];
        }
        
        // Reverse the SA12 and seperate mod1 and mod2
        if (SA12[i] % 3 == 1)
            s12[SA12[i]/3] = name;
        else 
            s12[SA12[i]/3 + n0] = name;
    }

    // Create Suffix Array for s12
	// recursion if names are not unique
    if (name == n02) {
        for (int i = 0; i < n02; i++)
            SA12[s12[i]-1] = i;
    }
    else {
        SA12.clear();
        skewAlgorithm(s12, SA12, n02, name);
        for (int i = 0; i < n02; i++)
            s12[SA12[i]] = i + 1;
    }

    // Construct suffix array for s0 and sort
    vector<uint32_t> SA0(n0);
    for (int i = 0, j = 0; i < n02; ++i)
        if (SA12[i] < n0)
            SA0[j++] = 3*SA12[i];
    radixSort(SA0, s, n0, 1, alphabet_size);

    // Merge SA0 and SA12 suffixes
    SA = vector<uint32_t>(n, 0);
    for (size_t p = 0, t = (n0-n1), k = 0; k < n; ++k) {

    const size_t i = SA12[t] < n0 ? SA12[t]*3 + 1 : (SA12[t] - n0)*3 + 2;
    const size_t j = SA0[p];

    const bool compare_pairs = SA12[t] < n0;
    const bool mod12_is_smaller =
      (compare_pairs ?
       leq_pair(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) :
       leq_triple(s[i], s[i+1], s12[SA12[t] + 1 - n0],
                  s[j], s[j+1], s12[j/3 + n0]));

    if (mod12_is_smaller) {
      SA[k] = i;
      ++t;
      if (t == n02)
        for (k++; p < n0; p++, k++)
          SA[k] = SA0[p];
    }
    else {
      SA[k] = j;
      ++p;
      if (p == n0)
        for (++k; t < n02; ++t, ++k)
          SA[k] = (SA12[t] < n0 ?
                   SA[k] = SA12[t]*3 + 1 :
                   (SA12[t] - n0)*3 + 2);
    }
  }
}

vector<int> kasai(vector<uint32_t> txt, vector<uint32_t> sa)
{
    int n = sa.size();
 
    vector<int> lcp(n, 0);
    vector<int> inv_sa(n, 0);
 
    for (int i=0; i < n; i++)
        inv_sa[sa[i]] = i;
 
    int k = 0;
    for (int i=0; i<n; i++)
    {
        if (inv_sa[i] == n-1)
        {
            k = 0;
            continue;
        }
 
        int j = sa[inv_sa[i]+1];
 
        while (i+k<n && j+k<n && txt[i+k]==txt[j+k])
            k++;
 
        lcp[inv_sa[i]] = k;
 
        if (k>0)
            k--;
    }
 
    return lcp;
}

int main(int argc, const char **argv){
    
    const string infile(argv[1]);
//    const string outfile(argv[2]);

    // Read file    
    string raw_s;
    readFasta(infile, raw_s);

    // Convert alphabet to number
    vector<uint32_t> s;
    for (auto i : raw_s) s.push_back(alphabetToNumerical(i));
    
    // Get starting timepoint
    auto start = std::chrono::high_resolution_clock::now();

    int initial_alphabet_size = 5;
    vector<uint32_t> SA;
    skewAlgorithm(s, SA, s.size(), initial_alphabet_size);

    vector<int> lcp;
    lcp = kasai(s, SA);

    // Get ending timepoint
    auto stop = std::chrono::high_resolution_clock::now();
 
    // Get duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    cout << endl;
    cout << "> Input sequence: ";
    for (char i : raw_s) cout << i; cout << endl;
    cout << "> Suffix array:" << endl;
    for (int i : SA) cout << std::setw(2) << i << " "; cout << endl;
    cout << "> LCP:" << endl;
    for (int i : lcp) cout << std::setw(2) << i << " "; cout << endl;

    cout << endl << "Execution time: " << duration.count() << " ms" << endl;
    
    return 0;
}