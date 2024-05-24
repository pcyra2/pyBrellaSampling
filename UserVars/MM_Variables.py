MM_DefaultVars = {  "parmfile" : "complex.parm7",
                    "ambercoor": "start.rst7",
                    "bincoordinates" : None,
                    "DCDfile" : "output",
                    "DCDfreq" : 1,
                    "restartname" : "output",
                    "restartfreq" : 1,              # Energy output and calculation also printed at this rate.
                    "outputname" : "output",
                    "outputTiming" : "100", 
                    "amber" : "on",
                    "switching" : "off",
                    "exclude" : "scaled1-4",
                    "1-4scaling" : 0.833333333,
                    "scnb" : 2.0,
                    "readexclusions" : "yes",
                    "cutoff" : 8.0,
                    "watermodel" : "tip3",
                    "pairListdist" : 11,
                    "LJcorrection" : "on",
                    "ZeroMomentum" : "off",
                    "rigidBonds" : "all",
                    "rigidTolerance" : "1.0e-8",
                    "rigidIterations" : 100,
                    "timeStep" : 2,                 # In fs
                    "fullElectFrequency" : 1,
                    "nonBondedFreq" : 1,
                    "stepspercycle" : 10,
                    "PME" : "off",
                    "PMEGridSizeX" : 300,
                    "PMEGridSizeY" : 300,
                    "PMEGridSizeZ" : 300,
                    "PMETolerance" : "1.0e-6",
                    "PMEInterpOrder" : 4,
                    "cellBasisVector1" : "135.913174 0.0 0.0",
                    "cellBasisVector2" : "-45.30439133333333 128.1401693173162 0.0",
                    "cellBasisVector3" : "-45.30439133333333 -64.0700846586581 -110.97264187403509",
                    "cellOrigin" : "0 0 0",
                    "cellBasisVector" : 135.913174, #### Base if all 3 are calculated...
                    "langevin" : "on",
                    "langevinDamping" : 5,
                    "langevinTemp" : 300,
                    "langevinHydrogen" : "off",
                    "temperature" : 300,
                    "BerendsenPressure" : "off",
                    "qmForces" : "off",
                    "CUDASOAintegrate" : "off",
                    "run" : 1000,
                    "minimize" : 0,
                    "qmParamPDB" : "syst-qm.pdb",
                    "qmColumn" : "beta",
                    "qmBondColumn" : "occ",
                    "QMsimsPerNode" : 1,
                    "QMElecEmbed" : "on",
                    "QMSwitching" : "on",
                    "QMSwitchingType" : "shift",
                    "QMPointChargeScheme" : "round",
                    "QMBondScheme" : "cs",
                    "qmBaseDir" : "/dev/shm/RUNDIR",
                    "qmConfigLine" : "! PBE 6-31G* EnGrad D3BJ TightSCF \n output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end",
                    "qmMult" : "1 1",
                    "qmCharge" : "1 0",
                    "qmSoftware" : "orca",
                    "qmExecPath" : "~/Software/ORCA/orca",
                    "QMOutStride" : 1,
                    "qmEnergyStride" : 1,
                    "QMPositionOutStride" : 1,
                    "colvars" : "off",
                    "colvarsConfig" : "colvars.conf",
                    }