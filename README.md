# Epilepsy-Jitter-Code

This repository contains code written for seizure spike detection and jitter analysis for Professor Yevgeny Berdichevsky's lab (https://www.lehigh.edu/~yeb211/) at Lehigh University. See the abstract of the work complete for this lab at the bottom of this ReadMe.

I recommend reading the Spike_Detection_User_Guide.docx before getting started.


# To run this code. 

I recommend reviewing the user guide. Quick start instructions are provided below.

1. Download all .m files into the same directory.

2. Open Seizure_Detection_And_Jitter_Analysis_Main in MATLAB. Run the code.

3. Select a .mat recording file.

4. Select spike time data (optional)

5. The data will be graphed onto the GUI. Take the following steps in order to record spike times of the seizure waveform.

    a. Zoom into one seizure waveform.

    b. Click "Find Spikes (algorithm)" button to find spikes in the waveform.

    c. If there are spikes the algorithm missed, click "Add spike" and click a point on the graph to add one additional spike time.

    d. If there are spikes that the algorithm incorrectly added, click "Remove Spike" and then click on the figure window the location of the spike to remove. 

    e. To set an alignment spike, click "Set Alignment Spike" button on the bottom right of the GUI. Click a point on the graph that most closely represents the alignment spike.

    f. If you zoom in or out on the figure, the spike time labels will be off. You can "Re-plot Peak Labels" to make the arrows match again.

    g. In order to save the spike times for the seizure, click the Save Burst button. It should then be added to the Bursts list. You can edit the burst again by double            clicking on it in the list. 

    h. Repeat the above steps to find more spike times of a different seizure.

    i. You can change the Find spikes algorithm parameters on the bottom if a seizure has lower/higher voltage than typical seizures.
    
    j. Save the spike time data to an excel file by filling out a filename and clicking Save Analysis.

    k. Close out of the window. The spike times and peak values are returned. Each column is a new seizure/burst analyzed and each row is a new spike.

6. A new window opens that allows you to find the jitter of the spike times recorded.

7. To find the jitter of different spike times than the ones just recorded. Click the button "Choose Spike Time Data"


# Description of what each of the .m files accomplishes:

Seizure_Detection_And_Jitter_Analysis_Main.m - main file user runs to perform spike detection and seizure jitter analysis.

selectspikes_gui.m: Graphical user interface to manually edit, add, and remove spikes. You can set an alignment spike (spike that aligns two seizure waveforms) as well as use the spike detection algorithm.

selectspike_gui.fig: The figure needed for the graphical user interface.

plot_peakpoints.m: Plots arrows at spike time points - used for visual inspection of the spike detection algorithm.

remove_peakpoints.m: Plots a white arrow over the black arrow to visually remove it on the figure.

findspikes.m: Uses the lower threshold, upper threshold, min spike period, max spike period, and min interspike interval to find spikes in the seizure waveform.

filter_data.m: Uses a band-stop filter to filter out 60 Hz, and a 4th order Butterworth high pass filter.

jitteranalysis_gui.m: Graphical user interface to find the jitter between spike times of different bursts.

jitteranlaysis_gui.fig: The figure needed for the jitter analysis graphical user interface.

MATLAB recording files: MEA_recording.mat, iMW_Recording.mat, sMW_Recording.mat

# Abstract:

Hippocampal and cortical slice-based models are widely used to study seizures and epilepsy. Seizure detection and quantification in these
models is an essential component for studies of mechanisms of epilepsy and for assessing therapeutic interventions. In order to obtain meaningful and
useful signals and maximize experimental throughput, variability should be minimized. Some electrical recording methods for detection of seizures
require insertion of an electrode into neuronal tissue, change in slice chemical microenvironment, and transients in temperature and pH. These
perturbations can cause acute and long-term alterations of the neuronal network which may be reflected in the variability of the recorded signal. In
this study we investigated the effect of experimental perturbations in three local field potential (LFP) recording methods including the substrate
micro-wires (s-MWs), multiple electrode arrays (MEAs), and inserted micro wire electrodes (i-MW). Use of these three methods enabled us to
isolate effects of different perturbations. We used organotypic hippocampal slices (OHCs) as an in-vitro model of posttraumatic epilepsy. To
investigate to what extent the paroxysmal events are affected by the disturbances caused by the recording method, we introduced jitter analysis,
which is sensitive to small differences in the seizure spike timing. Higher values of jitter indicate increased seizure variability. The use of jitter
analysis enabled us to track changes in variability on a relatively fine-grained time scale. We found that placement of a slice into new recording
medium can introduce long-lasting perturbations, while electrode insertion increased variability on a shorter time scale. OHCs also underwent state
transitions that occurred spontaneously, and that were characterized by transient increases in variability. The presence of these spontaneous transients
indicates that at least some of the variability is independent of the electrical recording method.

