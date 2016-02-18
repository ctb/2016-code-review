#include <iostream>
#include <string>
#include <cassert>

//
// these are givens - let's assume this code works.
//

typedef unsigned long long HashType;

#define uniqify_rc(f, r) ((f) < (r) ? (f) : (r))

#define twobit_repr(ch) ((ch) == 'A' ? 0LL : \
                         (ch) == 'T' ? 1LL : \
                         (ch) == 'C' ? 2LL : 3LL)

#define twobit_comp(ch) ((ch) == 'A' ? 1LL : \
                         (ch) == 'T' ? 0LL : \
                         (ch) == 'C' ? 3LL : 2LL)

HashType _hash(const char * kmer, const char k, HashType& _h, HashType& _r)
{
  HashType h = 0, r = 0;

  h |= twobit_repr(kmer[0]);
  r |= twobit_comp(kmer[k-1]);

  for (int i = 1, j = k - 2; i < k; i++, j--) {
    h = h << 2;
    r = r << 2;

    h |= twobit_repr(kmer[i]);
    r |= twobit_comp(kmer[j]);
  }

  _h = h;
  _r = r;

  return uniqify_rc(h, r);
}

void increment_count(HashType kmer)
{
  // stub code, for demo purposes
  std::cout << kmer << "\n";
}

class KmerIterator {
public:
  HashType bitmask;
  HashType forward, reverse;
  std::string _dna;
  int _ksize;
  int idx;
  
  KmerIterator(int ksize,
               std::string dna_string) {

    forward = 0;
    reverse = 0;
    bitmask = 0;
    _dna = dna_string;
    _ksize = ksize;
    idx = 0;
    
    // now fill bitmask with ones out to 2*k
    for (unsigned char i = 0; i < ksize; i++) {
      bitmask = (bitmask << 2) | 3;
    }
  }

  HashType first() {
    idx = _ksize;
    return _hash(_dna.c_str(), _ksize, forward, reverse);
  }

  HashType next() {
    char ch = _dna[idx];
    idx++;
    
    // take the last character in the window, which is sp[i];
    // add it to the hash of the forward strand.
    forward = forward << 2;
    forward |= twobit_repr(ch);
    forward &= bitmask;

    // add to the hash of the reverse strand.
    reverse = reverse >> 2;
    reverse |= (twobit_comp(ch) << (_ksize*2 - 2));

    // choose a unique representation of this k-mer
    return uniqify_rc(forward, reverse);
  }
  
};

//
// this is the function we're going to be looking at today.
// it's an optimized version of a function that extracts all
// DNA substrings of length k, calculates their forward and
// reverse-complement hashes, chooses the lower of the two, and
// passes that value to increment_count.
//

void count_kmers(const std::string &s, const unsigned int ksize)
{
  const char * sp = s.c_str();

  KmerIterator it(ksize, s);

  HashType hash = it.first();
  increment_count(hash);

  // start at end of first window, go to end of sequence.
  for (unsigned int i = ksize; i < s.length(); i++) {
    hash = it.next();
    increment_count(hash);
  }
}

int main()
{
  count_kmers("ATGGGACCAGATAGAGCCAGAGGACACATTAGGACGAT", 3);
}
