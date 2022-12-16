# MultibodyNonlinearBeams

##To run the Hamiltonian Neural Networks for the pendulum example:

1 - HNN/experiment_pend/data.py - where you can modify pendulum EoM or simply change the physical properties of the system if needed;

2 - HNN/experiment_pend/train.py - must run this code twice to generate baseline and hnn models, one setting the line 26 to action='store_false' 
and the other one to action='store_true';

3 - HNN/analyze-pend.ipynb - run the analysis in the Jupyter notebook.


##To run the D-HNN and MLP models associated with Matlab simulations for the beam FOM:

1 - SimulationFramework/beam_data_test.xlsx - where you can specify the properties of the beam;

1.1 - If you prefer you can define stiffness properties directly in lines 133-136 of SimulationFramework/background/create_beam.m;

2 - The structural damping coefficient of the model can be modified in line 11 of SimulationFramework/structure_properties.m;

3 - If you want to modify specific elements of the stiffness or damping matrices (e.g. to unclamp the beam and make it behave like a pendulum), 
that can be done uncommenting lines 45-50 as needed in MLP/data_maker.py;

4 - The ML models creation using the Matlab simulations, training and testing are done in MLP/MAIN.py. All settings for the Matlab simulations 
and training process can also be modified in this code. In line 173 you can save a .mat file with the results usign the different models.
