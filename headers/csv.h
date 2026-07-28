#ifndef CSV_H
#define CSV_H

#include <fstream>
#include <sstream>
#include <iostream>

#include <string>
#include <vector>

class CSV_Reader{
    public:
        std::string file_name;

        CSV_Reader(std::string file_name){
            this->file_name = file_name;
        }

        std::vector<std::vector<std::string>> getDataVector() {
            std::ifstream file(file_name);
            if (!file.is_open()) {
                std::cerr << "Couldn't read file: " << file_name << std::endl;
                return {};
            }

            std::vector<std::vector<std::string>> csvRows;
            std::string line;
            size_t lineNumber = 0;

            while (std::getline(file, line)) {
                ++lineNumber;
                if (lineNumber % 100 != 0) {
                    continue;
                }

                std::istringstream ss(line);
                std::vector<std::string> row;
                std::string value;

                while (std::getline(ss, value, ',')) {
                    row.push_back(value);
                }

                csvRows.push_back(std::move(row));
            }

            return csvRows;
        }
};

class CSV_Writer {
    public:
        std::string file_name;
        std::ofstream csvFile;

        CSV_Writer(std::string file_name){
            this->file_name = file_name;
            csvFile.open(file_name);
            if (!csvFile.is_open()) {
                std::cerr << "Error: Failed to open file '" << file_name << "' for writing." << std::endl;
             }   
        }

        void write_line(std::string data){
            csvFile << data + "\n";
        }

        void close_file(){
            csvFile.close();
        }
};
#endif