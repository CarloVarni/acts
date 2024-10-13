#pragma once

#include <unordered_set>
#include <utility>

namespace Acts {

  class DoubletCache {
  private:
    using key_t = std::pair<std::size_t, std::size_t>;
    struct HashFunction {
      std::size_t operator()(const key_t& key) const {
	return ( key.first << 32 ) | key.second;
      }
    };

    struct EqualFunction {
      bool operator()(const key_t& a, const key_t& b) const {
	return a.first == b.first and a.second == b.second;
      }
    };
    
    using cache_container_t = std::unordered_set< key_t,
						  Acts::DoubletCache::HashFunction,
						  Acts::DoubletCache::EqualFunction>;

  public:
    void insert(key_t&& a) {
      ++m_requests;
      auto [ptr, result] = m_cache.insert(std::move(a));
      if (not result) ++m_doppioni;
    }
    bool find(const key_t& a) const { return m_cache.find(a) != m_cache.end(); }    
    std::size_t size() const { return m_cache.size(); }

    std::size_t doppioni() const { return m_doppioni; }
    std::size_t requests() const { return m_requests; }

  private:
    cache_container_t m_cache{};
    std::size_t m_doppioni{};
    std::size_t m_requests{};
  };

}
