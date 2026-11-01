/*
RECONSTRUCCION basada en las salidas originales del usuario.
Compatible con PhysiCell 1.14.1; NO es un .cpp original recuperado.

Datos observados en los archivos XML/MAT originales:
  control: kcnj15_voltage=-30, loc401317_active=0;
  rosi:    kcnj15_voltage=-70, loc401317_active=0, salida G0/G1=0;
  CRISPRa: kcnj15_voltage=-70, loc401317_active=1,
           tasa de apoptosis guardada=9e9 por minuto.

El umbral -60 mV es HIPOTETICO: no se puede recuperar de estas trayectorias.
La implicacion causal voltaje->arresto/muerte es una regla asumida.
9e9/min reproduce un parametro historico numericamente saturado: fuerza muerte
en el motor; NO es una tasa biologica medida ni evidencia de eficacia.
No hay celulas inmunitarias ni eliminacion inmunitaria mecanistica.
La etiqueta loc401317_active no demuestra la identidad de CREB5-AS1.
No se deriva ningun parametro de gene_file.fa ni de results.tsv.
*/
#include "custom.h"
#include <stdexcept>

namespace {
int voltage_index=-1, activation_index=-1, apoptosis_index=-1;
bool use_rosiglitazone=false, use_crispra=false;
const double basal_voltage_mV=-30.0;       // observado en salidas del control
const double treated_voltage_mV=-70.0;     // observado en ambos brazos tratados
const double voltage_threshold_mV=-60.0;   // HIPOTESIS; umbral no identificado
const double historical_apoptosis_per_min=9.0e9; // observado; saturacion numerica

void check(bool ok,const std::string& message) {
    if (!ok) throw std::runtime_error("AML reconstruction: "+message);
}
}

void setup_microenvironment() { initialize_microenvironment(); }

void create_cell_types() {
    // v1.14.1 lee la semilla desde options/random_seed durante la carga XML.
    // No volver a sembrar aqui ni desde un callback paralelo.
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);
    cell_defaults.functions.volume_update_function=standard_volume_update_function;
    cell_defaults.functions.update_velocity=standard_update_cell_velocity;
    cell_defaults.functions.update_migration_bias=NULL;
    cell_defaults.functions.update_phenotype=NULL;
    cell_defaults.functions.custom_cell_rule=NULL;
    cell_defaults.functions.contact_function=NULL;
    cell_defaults.functions.add_cell_basement_membrane_interactions=NULL;
    cell_defaults.functions.calculate_distance_to_membrane=NULL;

    initialize_cell_definitions_from_pugixml();
    build_cell_definitions_maps();
    check(cell_definitions_by_index.size()==1,"Expected one cell definition");
    check(cell_defaults.phenotype.cycle.model().code==
        PhysiCell_constants::flow_cytometry_separated_cycle_model,
        "Use original separated flow cytometry cycle, code 6");
    voltage_index=cell_defaults.custom_data.find_variable_index("kcnj15_voltage");
    activation_index=cell_defaults.custom_data.find_variable_index("loc401317_active");
    check(voltage_index>=0 && activation_index>=0,"Missing original custom_data variables");
    apoptosis_index=cell_defaults.phenotype.death.find_death_model_index(
        PhysiCell_constants::apoptosis_death_model);
    check(apoptosis_index>=0,"Apoptosis death model missing");
    check(parameters.bools.find_index("terapia_rosiglitazona")>=0 &&
          parameters.bools.find_index("terapia_crispra")>=0,"Treatment flags missing");
    use_rosiglitazone=parameters.bools("terapia_rosiglitazona");
    use_crispra=parameters.bools("terapia_crispra");
    check(!(use_rosiglitazone && use_crispra),
        "Original data contain separate treatment arms, not a combination");
    setup_signal_behavior_dictionaries();
    setup_cell_rules();
    check(!PhysiCell_settings.rules_enabled,"Disable CBHG rules for this reconstruction");

    // Registrar DESPUES de leer las definiciones para que todas las celulas
    // hereden el callback correcto. No se modifica el nucleo de PhysiCell.
    cell_defaults.functions.update_phenotype=phenotype_function;
    display_cell_definitions(std::cout);
    if (use_crispra) {
        std::cout << "HISTORICAL REPLAY: apoptosis=9e9/min forces death; "
                     "not a calibrated therapeutic response.\n";
    }
}

void setup_tissue() {
    check(parameters.ints.find_index("number_of_cells")>=0,"number_of_cells missing");
    const int count=parameters.ints("number_of_cells");
    check(count>0,"number_of_cells must be positive");
    const std::vector<double>& box=microenvironment.mesh.bounding_box;
    double zmin=box[2], zmax=box[5];
    if (default_microenvironment_options.simulate_2D) zmin=zmax=0.0;
    // Rutina de siembra de la plantilla oficial. Las salidas originales parten
    // de cinco celulas. Mantiene incluso el sorteo de Z en el caso 2D.
    for (int n=0;n<count;++n) {
        std::vector<double> position(3,0.0);
        position[0]=box[0]+UniformRandom()*(box[3]-box[0]);
        position[1]=box[1]+UniformRandom()*(box[4]-box[1]);
        position[2]=zmin+UniformRandom()*(zmax-zmin);
        Cell* cell=create_cell(cell_defaults);
        cell->assign_position(position);
    }
    load_cells_from_pugixml();
    set_parameters_from_distributions();
    // No adelantar el callback: initial.xml conserva el valor XML original
    // -03.0, que significa -3.0, NO -30.0. Se documenta sin corregir el original.
}

void phenotype_function(Cell* cell,Phenotype& phenotype,double /* dt */) {
    // Correccion deliberada respecto al posible fallo historico: NO bloquear
    // la salida de la fase apoptotica escribiendo tasas de ciclo vivo en muertos.
    // Esto puede cambiar la persistencia de cadaveres del brazo rosi.
    if (phenotype.death.dead) return;

    double& voltage=cell->custom_data[voltage_index];
    double& active=cell->custom_data[activation_index];
    voltage=(use_rosiglitazone || use_crispra) ? treated_voltage_mV : basal_voltage_mV;
    active=use_crispra ? 1.0 : 0.0;

    // TRANSICION DE ESTADO: G0/G1 -> S en el ciclo nativo code=6.
    // Los otros enlaces siguen usando las duraciones originales del XML.
    // Esta condicion de voltaje reconstruye el comportamiento observado,
    // pero no identifica la condicion exacta que habia en el .cpp perdido.
    if (voltage<=voltage_threshold_mV && active<0.5) {
        phenotype.cycle.data.transition_rate(0,1)=0.0;
    }

    // TRANSICION A APOPTOSIS: requiere activacion y voltaje hiperpolarizado.
    // Se reproduce la tasa historica guardada. PhysiCell sortea muerte con
    // rate*dt: aqui supera 1 ampliamente y por tanto el evento es forzado.
    // El motor gestiona despues el ciclo apoptotico, volumen y eliminacion.
    if (voltage<=voltage_threshold_mV && active>0.5) {
        phenotype.death.rates[apoptosis_index]=historical_apoptosis_per_min;
    }
}

std::vector<std::string> my_coloring_function(Cell* cell) {
    return paint_by_number_cell_coloring(cell);
}
