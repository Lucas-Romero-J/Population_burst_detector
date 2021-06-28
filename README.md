# Population_burst_detector
Matlab script for population burst detection in MEA recordings from spontaneously active neurones

This code was developed in the context of microelectrode arrays (MEA) recordings in the superficial
dorsal horn of spinal cord. 
The paper reference will be included once published.

We recorded spontaneusly active neurones in an horizontal slice in vitro preparation, where we were
able to detect coordinated activity in which several units from the recordings participated.

With this script this population bursts can be detected and analised using a threshold criterium.
Other relevant data as the percentage of spikes falling within population bursts for each neuron,
the spike distribution in these bursts, the approximate duration or the bursts frequenzy can also
be extracted.

The user experience is quite simple. Once you run the script in Matlab there will be an file selection
window, and after visualizing the mean frequenzy distribution a horizontal threshold can be set, 
avoiding the background firing noise and allowing the detection of population bursts.

The relevant data for analysis purposes can be found in the generated output_file.csv and in the variables
percentages_units, Output_per_neurone and Experiment_summary

An example file containing the recording of 20 different neurons is provided for testing.
Supplementary files with raw data used in the publication are also included in the repository.
