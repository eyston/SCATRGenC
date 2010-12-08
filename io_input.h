#ifndef IO_INPUT_H_
#define IO_INPUT_H_

#include "structures.h"

SimulationParameters io_input_params(const char *input_file_name);

NBodies io_input_planets(const char *input_file_name);


#endif /* IO_INPUT_H_ */
