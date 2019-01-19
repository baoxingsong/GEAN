//
// Created by Baoxing song on 10.08.18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include <bitset>
TEST(bitset, c1){
    std::bitset<13> bitTest;
    bitTest[0] = 0;
    bitTest[0] = 1;
    bitTest[1] = 1;
    std::cout << std::endl << "foo: " << bitTest << std::endl;
}

TEST(bool, c1){
    bool b = -1;
    std::cout << std::endl << "b: " << b << std::endl;
    std::cout << std::endl << "0-b: " << 0-b << std::endl;
    b = false;
    std::cout << std::endl << "b: " << b << std::endl;
    std::cout << std::endl << "0-b: " << 0-b << std::endl;
}
