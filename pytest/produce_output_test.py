import pytest

import os
import subprocess

common_init_params = """
/PhysicsList/RegisterPhysics G4EmStandardPhysics_option4
/PhysicsList/RegisterPhysics G4DecayPhysics
/PhysicsList/RegisterPhysics G4RadioactiveDecayPhysics
/PhysicsList/RegisterPhysics G4OpticalPhysics
/PhysicsList/RegisterPhysics NexusPhysics
/PhysicsList/RegisterPhysics G4StepLimiterPhysics

/nexus/RegisterTrackingAction DefaultTrackingAction
/nexus/RegisterEventAction DefaultEventAction
/nexus/RegisterRunAction DefaultRunAction

/nexus/RegisterPersistencyManager PersistencyManager

"""

single_part_params = """
/Generator/SingleParticle/particle e-
/Generator/SingleParticle/min_energy 10. keV
/Generator/SingleParticle/max_energy 10. keV

"""

next100_params = """
/Geometry/Next100/elfield true
/Geometry/Next100/EL_field 13 kV/cm
/Geometry/Next100/max_step_size 1. mm
/Geometry/Next100/pressure 15. bar
/Geometry/Next100/sc_yield 10000 1/MeV

"""

def run_simulation(NEXUSDIR, init_path):
    my_env    = os.environ
    nexus_exe = NEXUSDIR + '/bin/nexus'
    command   = [nexus_exe, '-b', '-n', '1', init_path]
    p         = subprocess.run(command, check=True, env=my_env)

@pytest.mark.order(1)
def test_create_nexus_output_file_next100(config_tmpdir, output_tmpdir,
                                          NEXUSDIR,
                                          full_base_name_next100,
                                          nexus_full_output_file_next100):
    # Init file
    init_text = f"""
/nexus/RegisterGeometry Next100OpticalGeometry

/nexus/RegisterGenerator SingleParticleGenerator

/nexus/RegisterMacro {config_tmpdir}/{full_base_name_next100}.config.mac
"""

    init_text = f'{common_init_params} {init_text}'
    init_path = os.path.join(config_tmpdir, full_base_name_next100+'.init.mac')
    init_file = open(init_path,'w')
    init_file.write(init_text)
    init_file.close()

    # Config file
    config_text = f"""
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/process/em/verbose 0

/Generator/SingleParticle/region CENTER

/nexus/persistency/save_strings true
/nexus/persistency/output_file {output_tmpdir}/{full_base_name_next100}
/nexus/random_seed 21051817
"""
    config_text = f'{config_text} {next100_params} {single_part_params}'
    config_path = os.path.join(config_tmpdir, full_base_name_next100+'.config.mac')
    config_file = open(config_path,'w')
    config_file.write(config_text)
    config_file.close()

    # Running the simulation
    run_simulation(NEXUSDIR, init_path)

    return nexus_full_output_file_next100



@pytest.mark.order(2)
def test_create_nexus_output_file_new(config_tmpdir, output_tmpdir, NEXUSDIR,
                                      full_base_name_new,
                                      nexus_full_output_file_new):
    # Init file
    init_text = f"""
/nexus/RegisterGeometry NextNew

/nexus/RegisterGenerator SingleParticleGenerator

/nexus/RegisterMacro {config_tmpdir}/{full_base_name_new}.config.mac
"""
    init_text = f'{common_init_params} {init_text}'
    init_path = os.path.join(config_tmpdir, full_base_name_new+'.init.mac')
    init_file = open(init_path,'w')
    init_file.write(init_text)
    init_file.close()

    # Config file
    config_text = f"""
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/process/em/verbose 0

/Geometry/NextNew/elfield true
/Geometry/NextNew/EL_field 13 kV/cm
/Geometry/NextNew/max_step_size 1. mm
/Geometry/NextNew/pressure 15. bar
/Geometry/NextNew/sc_yield 10000 1/MeV

/Generator/SingleParticle/region CENTER

/nexus/persistency/save_strings true
/nexus/persistency/output_file {output_tmpdir}/{full_base_name_new}
/nexus/random_seed 21051817
"""
    config_text = f'{config_text} {single_part_params}'
    config_path = os.path.join(config_tmpdir, full_base_name_new+'.config.mac')
    config_file = open(config_path,'w')
    config_file.write(config_text)
    config_file.close()

    # Running the simulation
    run_simulation(NEXUSDIR, init_path)

    return nexus_full_output_file_new



@pytest.mark.order(3)
def test_create_nexus_output_file_flex100(config_tmpdir, output_tmpdir, NEXUSDIR,
                                          full_base_name_flex100,
                                          nexus_full_output_file_flex100):
    # Init file
    init_text = f"""
/nexus/RegisterGeometry NextFlex

/nexus/RegisterGenerator SingleParticleGenerator

/nexus/RegisterMacro {config_tmpdir}/{full_base_name_flex100}.config.mac
"""
    init_text = f'{common_init_params} {init_text}'
    init_path = os.path.join(config_tmpdir, full_base_name_flex100+'.init.mac')
    init_file = open(init_path,'w')
    init_file.write(init_text)
    init_file.close()

    # Config file
    config_text = f"""
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/process/em/verbose 0

/Geometry/NextFlex/gas                   enrichedXe
/Geometry/NextFlex/gas_pressure          15. bar
/Geometry/NextFlex/gas_temperature       300. kelvin
/Geometry/NextFlex/sc_yield              25510. 1/MeV
/Geometry/NextFlex/e_lifetime            12. ms
/Geometry/NextFlex/active_length         1204.95 mm
/Geometry/NextFlex/active_diam           984.0 mm
/Geometry/NextFlex/drift_transv_diff     1. mm/sqrt(cm)
/Geometry/NextFlex/drift_long_diff       .2 mm/sqrt(cm)
/Geometry/NextFlex/buffer_length         254.6 mm
/Geometry/NextFlex/cathode_transparency  .98
/Geometry/NextFlex/anode_transparency    .88
/Geometry/NextFlex/gate_transparency     .88
/Geometry/NextFlex/el_gap_length         10.  mm
/Geometry/NextFlex/el_field_on           true
/Geometry/NextFlex/el_field_int          16. kilovolt/cm
/Geometry/NextFlex/el_transv_diff        0. mm/sqrt(cm)
/Geometry/NextFlex/el_long_diff          0. mm/sqrt(cm)
/Geometry/NextFlex/fc_wls_mat            TPB
/Geometry/NextFlex/fc_with_fibers        false
/Geometry/NextFlex/fiber_mat             EJ280
/Geometry/NextFlex/fiber_claddings       2
/Geometry/NextFlex/fiber_sensor_time_binning  25. ns
/Geometry/NextFlex/ep_with_PMTs          true
/Geometry/NextFlex/ep_with_teflon        false
/Geometry/NextFlex/ep_copper_thickness   12. cm
/Geometry/NextFlex/ep_wls_mat            TPB
/Geometry/PmtR11410/time_binning         25. ns
/Geometry/NextFlex/tp_copper_thickness   12. cm
/Geometry/NextFlex/tp_teflon_thickness   2.1 mm
/Geometry/NextFlex/tp_teflon_hole_diam   7. mm
/Geometry/NextFlex/tp_wls_mat            TPB
/Geometry/NextFlex/tp_kapton_anode_dist  15. mm
/Geometry/NextFlex/tp_sipm_sizeX         1.3 mm
/Geometry/NextFlex/tp_sipm_sizeY         1.3 mm
/Geometry/NextFlex/tp_sipm_sizeZ         2.0 mm
/Geometry/NextFlex/tp_sipm_pitchX        15.55 mm
/Geometry/NextFlex/tp_sipm_pitchY        15.55 mm
/Geometry/NextFlex/tp_sipm_time_binning  1. microsecond
/Geometry/NextFlex/ics_thickness         12. cm

/Generator/SingleParticle/region         AD_HOC
/Geometry/NextFlex/specific_vertex       0. 0. 500. mm

/nexus/persistency/save_strings true
/nexus/persistency/output_file {output_tmpdir}/{full_base_name_flex100}
/nexus/random_seed 21051817
"""
    config_text = f'{config_text} {single_part_params}'
    config_path = os.path.join(config_tmpdir, full_base_name_flex100+'.config.mac')
    config_file = open(config_path,'w')
    config_file.write(config_text)
    config_file.close()

    # Running the simulation
    run_simulation(NEXUSDIR, init_path)

    return nexus_full_output_file_flex100


@pytest.mark.order(4)
def test_create_nexus_output_file_demopp(config_tmpdir, output_tmpdir, NEXUSDIR, full_base_name_demopp, nexus_full_output_file_demopp):

    for run in ["run5", "run7", "run8", "run9", "run10"]:

        init_text = f"""
    /nexus/RegisterGeometry NextDemo

    /nexus/RegisterGenerator SingleParticleGenerator

    /nexus/RegisterMacro {config_tmpdir}/{full_base_name_demopp.format(run=run)}.config.mac
    """
        init_text = f'{common_init_params} {init_text}'
        init_path = os.path.join(config_tmpdir, full_base_name_demopp.format(run=run) +'.init.mac')
        init_file = open(init_path,'w')
        init_file.write(init_text)
        init_file.close()

        config_text = f"""
    /run/verbose 1
    /event/verbose 0
    /tracking/verbose 0

    /Geometry/NextDemo/config {run}
    /Geometry/NextDemo/elfield true
    /Geometry/NextDemo/EL_field_intensity 13 kV/cm
    /Geometry/NextDemo/max_step_size 1. mm
    /Geometry/NextDemo/pressure 10. bar
    /Geometry/NextDemo/sc_yield 10000 1/MeV

    /Geometry/NextDemo/specific_vertex 0. 0. 10. cm

    /Generator/SingleParticle/region AD_HOC

    /nexus/persistency/save_strings true
    /nexus/persistency/output_file {output_tmpdir}/{full_base_name_demopp.format(run=run)}
    /nexus/random_seed 21051817

    """
        config_text = f'{config_text} {single_part_params}'
        config_path = os.path.join(config_tmpdir, full_base_name_demopp.format(run=run) +'.config.mac')
        config_file = open(config_path,'w')
        config_file.write(config_text)
        config_file.close()

        # Running the simulation
        run_simulation(NEXUSDIR, init_path)

    return nexus_full_output_file_demopp


@pytest.mark.order(5)
def test_create_nexus_output_file_no_strings(config_tmpdir, output_tmpdir,
                                             NEXUSDIR,
                                             base_name_no_strings,
                                             nexus_output_file_no_strings):
    # Init file
    init_text = f"""
/nexus/RegisterGeometry Next100OpticalGeometry

/nexus/RegisterGenerator SingleParticleGenerator

/nexus/RegisterMacro {config_tmpdir}/{base_name_no_strings}.config.mac
"""
    init_text = f'{common_init_params} {init_text}'
    init_path = os.path.join(config_tmpdir, base_name_no_strings+'.init.mac')
    init_file = open(init_path,'w')
    init_file.write(init_text)
    init_file.close()

    # Config file
    config_text = f"""
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/process/em/verbose 0

/Generator/SingleParticle/region CENTER

/nexus/persistency/save_strings false
/nexus/persistency/output_file {output_tmpdir}/{base_name_no_strings}
/nexus/random_seed 21051817
"""
    config_text = f'{config_text} {next100_params} {single_part_params}'
    config_path = os.path.join(config_tmpdir, base_name_no_strings+'.config.mac')
    config_file = open(config_path,'w')
    config_file.write(config_text)
    config_file.close()

    # Running the simulation
    run_simulation(NEXUSDIR, init_path)

    return nexus_output_file_no_strings