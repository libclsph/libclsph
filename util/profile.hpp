#ifndef _PROFILE_HPP_
#define _PROFILE_HPP_

#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>

#define profile(name, fun) if(profile_block p = profile_block(name, fun))

class profile_block {
public:

    //Don't assign default values
    enum steps_enum {
        STEP_1,
        STEP_2,
        SORT,
        MEMORY_TRANSFERS,
        FILE_IO
    };

    profile_block(std::string name, std::function<void(std::string, long long)> f = COUT_LOG)
        : block_name(name),
          func(f) {
        start = std::chrono::high_resolution_clock::now();
    }

    profile_block(steps_enum step_id, std::function<void(steps_enum, long long)> f = TALLY_STEP_TIME)
        : current_step(step_id),
          func_tally(f) {
        start = std::chrono::high_resolution_clock::now();
    }

    ~profile_block() {
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds =
            std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        if(func) {
            func(block_name, microseconds);
        }
        if(func_tally) {
            func_tally(current_step,microseconds);
        }
    }

    static void print_stats() {

        long long total_time = 0;

        for ( int step_ID = STEP_1; step_ID <= FILE_IO; step_ID++ ) {
            steps_enum step = static_cast<steps_enum>(step_ID);
            total_time +=step_run_length[step];
        }

        for ( int step_ID = STEP_1; step_ID <= FILE_IO; step_ID++ ) {
            steps_enum step = static_cast<steps_enum>(step_ID);

            std::cout << std::fixed << std::setprecision(3) <<
                      ((double)step_run_length[step]) / 1000.f << " milliseconds were spent on " <<
                      step_names[step] << " and accounted for " << ((double)step_run_length[step]/(double)total_time)*100 << "%" << " of the processing time." << std::endl;
        }

    }

    static void reset_stats(){
        profile_block::step_run_length = { { STEP_1, 0 },{ STEP_2, 0 },{ SORT, 0 },{ MEMORY_TRANSFERS, 0 },{ FILE_IO, 0 } };
    }

    operator bool() const {
        return true;
    }

    static std::function<void(std::string, long long)> COUT_LOG;
    static std::function<void(steps_enum, long long)> TALLY_STEP_TIME;

    static std::map<steps_enum,std::string> step_names;
    static std::map<steps_enum,long long> step_run_length;

private:
    std::chrono::high_resolution_clock::time_point start;
    std::string block_name;
    steps_enum  current_step;
    std::function<void(std::string, long long)> func;
    std::function<void(steps_enum, long long)> func_tally;
};

std::map<profile_block::steps_enum,std::string> profile_block::step_names    = { { STEP_1,           "Sph step 1      " },
                                                                                 { STEP_2,           "Sph step 2      " },
                                                                                 { SORT,             "Sort            " },
                                                                                 { MEMORY_TRANSFERS, "Memory transfers" },
                                                                                 { FILE_IO,          "File IO         " }};

std::map<profile_block::steps_enum,long long>   profile_block::step_run_length = { { STEP_1, 0 },{ STEP_2, 0 },{ SORT, 0 },{ MEMORY_TRANSFERS, 0 },{ FILE_IO, 0 } };

std::function<void(profile_block::steps_enum,long long)> profile_block::TALLY_STEP_TIME=
[] (profile_block::steps_enum step, long long elapsed) {
    step_run_length[step]+=elapsed;
};

std::function<void(std::string, long long)> profile_block::COUT_LOG =
[] (std::string block_name, long long elapsed) {
    std::cout << std::fixed << std::setprecision(3) <<
              "[PROFILE] " << block_name	<<
              " completed in " << ((double)elapsed) / 1000.f <<
              " milliseconds" << std::endl;
};

#endif
