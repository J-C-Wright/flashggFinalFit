



void pochtest() {

    double x = -1.99;
    unsigned n = 4;

    double result(x);
    for (unsigned i(1);i<n;i++) {
        std::cout << result << std::endl;
        result *= x + i;
    }

    std::cout << result << std::endl;

}
