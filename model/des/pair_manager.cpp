// example

#include <unordered_map>
#include <utility>
#include <iostream>

std::unordered_map<uint64_t,char> map;

extern "C" {
  void ht_add_pair(int ii, int jj) {
    uint64_t key = ((uint64_t(ii) << 32) | jj);
    // default value is zero, what we want
    map[key] = map[key] + 1;
    //    std::cout << " ADD~~~~" << ii << "   " << jj << " key: " << std::hex << key << " count is==" << int(map[key]) << "  " << sizeof(ii) << "  "<< sizeof(key) << std::endl;
    if (map[key] > 10) {
      std::cout << " ERROR --- A PARTICLE SHOULD NOT BE IN MORE THAN EIGHT SAPS " << map[key];
      std::exit(-1);
    }
  }

  void ht_del_pair(int ii, int jj) {
    uint64_t key = ((uint64_t(ii) << 32) | jj);
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
    uint64_t key = ((uint64_t(ii) << 32) | jj);
    // default value is zero, what we want
    std::cout << ( map.find(key) != map.end() ) << " CHECKING PAIR~~~~~~~~~~~~~~~~" << ii << "   " << jj << std::endl;
    return ( map.find(key) != map.end() );
  }

}
