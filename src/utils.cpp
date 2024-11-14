#include <cstdint>
#include <sstream>
#include <fstream>
#include "../include/utils.hpp"

namespace kaminari::utils {

std::vector<std::string>
read_filenames(std::string const& filenames_list) 
{
    std::vector<std::string> buffer;
    std::ifstream in(filenames_list);
    if (!in.is_open()) throw std::runtime_error("error in opening file");
    std::string filename;
    while (in >> filename) buffer.push_back(filename);
    in.close();
    return buffer;
}

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier)
{
    std::stringstream filename;
    filename 
    << tmp_dirname 
    << (tmp_dirname != "" ? "/" : "") 
    << prefix << "_" << run_identifier << ".bin";
    return filename.str();
}

uint64_t getTotalVirtualMem() {
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    uint64_t totalVirtualMem = memInfo.totalram;
    //Add other values in next statement to avoid int overflow on right hand side...
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;
    return totalVirtualMem/1000/1000;
}

uint64_t getVirtualMemUsed() {
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    uint64_t virtualMemUsed = memInfo.totalram - memInfo.freeram;
    //Add other values in next statement to avoid int overflow on right hand side...
    virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
    virtualMemUsed *= memInfo.mem_unit;
    return virtualMemUsed/1000/1000;
}

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}
int getVirtualMemUsedByProcess(){ //Note: this value is in MB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result/1000;
}

uint64_t getTotalPhysMem() {
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    uint64_t totalPhysMem = memInfo.totalram;
    //Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;
    return totalPhysMem/1000/1000;
}

uint64_t getPhysMemUsed() {
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    uint64_t physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;
    return physMemUsed/1000/1000;
}

int getPhysMemUsedByProcess(){ //Note: this value is in MB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result/1000;
}

void printRAMInfo() {
    std::cerr << "VirtualMem:\n\tTotal: " << getTotalVirtualMem() << " MB\tUsed: " << getVirtualMemUsed() << " MB\n\tUsedByProcess: " << getVirtualMemUsedByProcess() << " MB\n";
    std::cerr << "PhysMem:\n\tTotal: " << getTotalPhysMem() << " MB\tUsed: " << getPhysMemUsed() << " MB\n\tUsedByProcess (~RAM usage): " << getPhysMemUsedByProcess() << " MB\n";

}

} // namespace util