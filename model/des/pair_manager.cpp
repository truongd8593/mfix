// example

#include <unordered_map>
#include <utility>
#include <iostream>

std::unordered_map<long,char> map;

extern "C" {
  void ht_add_pair(int ii, int jj) {
    long key = ((ii << 32) & jj);
    // default value is zero, what we want
    map[key] = map[key] + 1;
    std::cout << " ADDED PAIR~~~~~~~~~~~~~~~~" << ii << "   " << jj << std::endl;
  }

  void ht_del_pair(int ii, int jj) {
    long key = ((ii << 32) & jj);
    // default value is zero, what we want
    if ( map.find(key) != map.end() ){
    map[key] = map[key] - 1;
    if (map[key == 0]) {
      map.erase(key);
      std::cout << " ERASED PAIR~~~~~~~~~~~~~~~~" << ii << "   " << jj << std::endl;
    }
    }
  }

  bool ht_is_pair(int ii, int jj) {
    long key = ((ii << 32) & jj);
    // default value is zero, what we want
    std::cout << ( map.find(key) != map.end() ) << " CHECKING PAIR~~~~~~~~~~~~~~~~" << ii << "   " << jj << std::endl;
    return ( map.find(key) != map.end() );
  }

}
