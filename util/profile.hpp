#ifndef _PROFILE_HPP_
#define _PROFILE_HPP_

#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>

class profile_data {
public:
    typedef size_t key_t;
    std::map<key_t, unsigned long long> data;
};

#define profile(name, fun) if(profile_block p = profile_block(name, fun))

class profile_block {
public:

    profile_block(profile_data::key_t bucket_key, profile_data* data)
        : profile_bucket_key(bucket_key),
        data(data) {

        start = std::chrono::high_resolution_clock::now();
    }

    ~profile_block() {
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds =
            std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

        data->data[profile_bucket_key] += microseconds;
    }

    operator bool() const {
        return true;
    }

private:
    std::chrono::high_resolution_clock::time_point start;
    profile_data::key_t profile_bucket_key;
    profile_data* data;
};

#endif
