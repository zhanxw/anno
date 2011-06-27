#include "LogFile.h"
#include <ctime>

int main(int argc, char *argv[])
{
#if 0
    LogFile::getLogger()->open(argv[0], "test.log");
    LogFile::getLogger()->close();
#else
    LOG_START("test.log");
    LOG_START_TIME;
    LOG << "abc";
    LOG << "abc"<<std::endl;;
    LOG_END_TIME;
    LOG_END ;
#endif
    return 0;
}
