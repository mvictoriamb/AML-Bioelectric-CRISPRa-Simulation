#ifndef AML_RECOVERED_CUSTOM_H
#define AML_RECOVERED_CUSTOM_H
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
using namespace BioFVM;
using namespace PhysiCell;
void create_cell_types();
void setup_microenvironment();
void setup_tissue();
void phenotype_function(Cell*, Phenotype&, double);
std::vector<std::string> my_coloring_function(Cell*);
#endif
