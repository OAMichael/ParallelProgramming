#ifndef COMMON_H
#define COMMON_H

#include <fstream>


#define ISIZE 5000
#define JSIZE 5000

template<typename T>
struct My2DArray {
    T* data;

    size_t xSize;
    size_t ySize;

    My2DArray(const size_t xS, const size_t yS) {
        data = new T[xS * yS];
        xSize = xS;
        ySize = yS;
    }

    ~My2DArray() { delete[] data; }

    inline T* operator[](const size_t idx) { return data + ySize * idx; }
    inline const T* operator[](const size_t idx) const { return data + ySize * idx; }
};


template<typename T>
void serializeToFile(const std::string& outFilename, const My2DArray<T>& a) {
    std::ofstream outFile;
    outFile.open(outFilename);
    if (!outFile.is_open()) {
        std::cerr << "Could not open file: " << outFilename << std::endl;
    }
    else {
        for(int i = 0; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                outFile << a[i][j] << " ";
            }
            outFile << std::endl;
        }
        outFile.close();
    }
}

#endif