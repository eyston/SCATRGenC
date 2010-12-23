#ifndef IO_INPUT_H_
#define IO_INPUT_H_

#include "structures.h"

SimulationParameters io_input_params(const char *input_file_name);

void io_input_planets(const char *input_file_name, size_t &nbod, size_t &npl, float mass[NPLMAX], float rpl[NPLMAX], float xh[NPLMAX], float yh[NPLMAX], float zh[NPLMAX], float vxh[NPLMAX], float vyh[NPLMAX], float vzh[NPLMAX]);

void io_input_particles(const char *input_file_name, size_t &ntp, float xht[NTPMAX], float yht[NTPMAX], float zht[NTPMAX], float vxht[NTPMAX], float vyht[NTPMAX], float vzht[NTPMAX]);

#endif /* IO_INPUT_H_ */
