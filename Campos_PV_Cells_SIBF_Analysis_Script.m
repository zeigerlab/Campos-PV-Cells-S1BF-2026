%% Data analysis from "Dynamic recruitment of PV interneurons in the somatosensory cortex induced by experience dependent plasticity"
% PMID: 40667001  PMCID: PMC12262626
% See manuscript for full methods
% Data analysis was performedin MATLAB_R2020b

% The initial sections of code are required for formatting and organization
% of data to generate subsequent plots. Any sections with "Required" in the
% section title must be run for subsequent plots to work. After the
% required sections of code are run, subsequent plots can be run
% independently as desired.
%% Section 1 (Required): Load the data tables
% This loads data tables for longitudinal imaging of PV interneurons from
% S1BF. PV cells expressed GCaMP6s and fluoresence traces over time were
% extracted using DC Calcium. Modified Z-scores were calculated and
% bootstrapping was performed to determine responsivity to passive whisker
% deflections. Additional metrics were quantified, as described below. 

clear %clear workspace variables

%Load the whisker evoked data table
% load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Iris/Chronic Whisker Trimming Experiment/Analysis - Will/DataByCell_wCorr.mat')
load('DataByCell_wCorr.mat')

% Each row of the table represents data from one PV cell in one imaging
% session. Within a given mouse, multiple areas (barrels) were imaged at
% two depths (Z-planes). Cells are assigned an index so that the same cell
% within a given mouse can be tracked across multiple imaging sessions.

%'MouseID' - ID of mouse

%'Timepoint' - Imaging timepoint relative to whisker trimming

%'Group' - Experiments were conducted on small groups of mice all imaged on
%the same dates

%'Whisker_stim' - Whisker deflected during imaging

%'Barrel_imaged' - Barrel over which the imaging FOV was centered

%'Gender' - Mouse sex (male/female)

%'Index' - Cell index within a mouse for longitudinal tracking.

%'Z_plane_Index' - Z plane in L2/3 imaged. 1=superficial, 2= deeper.

%'signal_trace' - Modified Z-scored fluoresence trace.

%'stim_epochs' - Time epochs during whisker stimulation set to 1.

%'activeROI' - whether the cell was active (1) or not (0) during imaging,
%defined as Z-score >3 for at least 0.5 s consecutively.

%'signal_percentile' - statistic from bootstrapping to determine if
%responses are time-locked to whisker stimulation. <0.01 considered
%responsive. The calcium trace is shuffled and the correlation coefficient
%between the real (or shuffled) trace is compared to the 'stim_epochs'
%trace.

%'MI_signal_percentile'

%'Responsivity' - set to 0 (not responsive) or 1 (responsive) based on
%'signal percentile'.

%'allstimsallROIs' - a matrix representing 20 individual whisker deflection
%epochs, 1 second before and ~3 seconds after onset of stimulation.

%'allstimsROIs_mean' - the mean trace of 'allstimsallROIs'

%'allstimsallROIs_AUC - the AUC of allstimsROIs_mean

%'PeaksMax_mean' - the mean amplitude of whisker evoked calcium transients,
%including only stimulations with a peak above criterion (exlcuding
%stimulations without an evoked peak

%'PeaksMax_meanallstims' - the mean amplitude of whisker evoked calcium transients,
%including all 20 stimulations even those stimulations without an evoked peak

%'PeaksLatency_mean' - the mean latency to peak of whisker evoked calcium
%transients

%'PrctStimActive' - the percentage of 20 whisker stimulations with an
%evoked peak transient

%'signal_percentile_stim' - statistic from bootstrapping to determine if
%responses are time-locked to whisker stimulation. <0.01 considered
%responsive. The 'stim_epochs' trace is shuffled and the correlation coefficient
%between the real (or shuffled) trace is compared to the 'signal_trace'

%'MI_signal_percentile_stim'

%'ROI_centers' - x and y position of the cell within the imaging FOV

%'MI' - modulation index of the cell

%'ResponsivityStim' - set to 0 (not responsive) or 1 (responsive) based on
%'signal_percentile_stim'.

%'MI_ResponsivityStim'

%Save the data with alternate names so subsequent analyses work properly
Data=Data_ByCell_wCorr;
Data_ByCell=Data_ByCell_wCorr;

%Remove the MI >1 threshold for responsivity so users can choose their own
%threshold or set back to 1
for i=1:size(Data,1)
    if Data.signal_percentile_stim{i} <= 0.01 && Data.activeROI(i)==1 %Responsive based on signal correlation with stimulus trace
        Data.ResponsivityStim(i) = 1;
    else
        Data.ResponsivityStim(i) = 0;
    end
end
clear i

%Load the no-trimming control data
% load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/Control_ExpDepPlast/Imaging Control Trimming/Analysis/PV_ExpDepPlas_NoTrim_DataByCell.mat')
load('PV_ExpDepPlas_NoTrim_DataByCell.mat')
Data_NoTrim=Data_ByCell_NoTrim;

%Load the spontaneous (non whisker evoked) data
% load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/PVG-EXPERIMENTS/2P Plasticity/Calcium Imaging/Spontaneous Data Analysis/Spont_DataByCell_wCorr.mat')
load('Spont_DataByCell_wCorr.mat')
Spont_Data_ByCell=Spont_Data_ByCell_wCorr;


%Trim evoked data to just one whisker to match one condition (spontaneous) for spontaneous data
data_core=Data_ByCell(Data_ByCell.Whisker_stim=='Stm_D1',:); 
%sort rows so spontaneous and evoked data are aligned
Spont_Data_ByCell=sortrows(Spont_Data_ByCell,{'MouseID','Timepoint','Barrel_imaged','Zplane'});
data_core=sortrows(data_core,{'MouseID','Timepoint','Barrel_imaged','Z_plane_Index'});
%Match the cell index from evoked data to spontaneous data
Spont_Data_ByCell.index=data_core.Index;

%% Fig S1b-d (Optional): Determine a MI threshold for classifying responsive cells
% This code will plot responsive cells with different MI thresholds (<0.5,
% 0.5-1, and 1-2) as in Fig. S1b-d.

C1_BL = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL",:);
C1_BL_resp = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 1,:);
C1_BL_non_resp = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 0,:);
D1_BL = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL",:);

% Get means for each barrel/stim/timepoint from individual mice using shuffled stimulus responsivity criteria
omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);
T1=Data(:,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim','Gender'});
[G,mean_ByMouse]=findgroups(T1);
C=groupcounts(T1,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim','Gender'});

mean_ByMouse.NumRespCells=splitapply(osum,Data.ResponsivityStim,G);
mean_ByMouse.TotalCells=C.GroupCount;
mean_ByMouse.PercentResponsive=(mean_ByMouse.NumRespCells./mean_ByMouse.TotalCells)*100;

mean_ByMouse.MI_NumRespCells=splitapply(osum,Data.MI_ResponsivityStim,G);
mean_ByMouse.MI_PercentResponsive=(mean_ByMouse.MI_NumRespCells./mean_ByMouse.TotalCells)*100;
%Then summarize by group to get mean with n=# of mice
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_ByGroup]=findgroups(T1);

mean_ByGroup.PercentResponsive=splitapply(omean,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.PercentResponsive_SEM=splitapply(sem,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.MI_PercentResponsive=splitapply(omean,mean_ByMouse.MI_PercentResponsive,G);
mean_ByGroup.MI_PercentResponsive_SEM=splitapply(sem,mean_ByMouse.MI_PercentResponsive,G);

% Pool all cells by group for responsive cells and non-responsive cells
Data_Resp=Data(Data.ResponsivityStim==1,:);
T1=Data_Resp(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_Resp]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
% sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,1,'omitnan')/sqrt(size(x,1)));

mean_Resp.TotalCells=C.GroupCount;
mean_Resp.AUC=splitapply(omean,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.AUC_SEM=splitapply(sem,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.MI=splitapply(omean,Data_Resp.MI,G);
mean_Resp.MI_SEM=splitapply(sem,Data_Resp.MI,G);
mean_Resp.Amp=splitapply(omean,Data_Resp.PeaksMax_mean,G);
mean_Resp.Amp_SEM=splitapply(sem,Data_Resp.PeaksMax_mean,G);
mean_Resp.Latency=splitapply(omean,Data_Resp.PeaksLatency_mean,G);
mean_Resp.Latency_SEM=splitapply(sem,Data_Resp.PeaksLatency_mean,G);
mean_Resp.PrctStimActive=splitapply(omean,Data_Resp.PrctStimActive,G);
mean_Resp.PrctStimActive_SEM=splitapply(sem,Data_Resp.PrctStimActive,G);
mean_Resp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp.allstimsROIs_mean),G);
mean_Resp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp.allstimsROIs_mean),G);

Data_NResp=Data(Data.ResponsivityStim==0,:);
T1=Data_NResp(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_NResp]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});

mean_NResp.TotalCells=C.GroupCount;
mean_NResp.AUC=splitapply(omean,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.AUC_SEM=splitapply(sem,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.MI=splitapply(omean,Data_NResp.MI,G);
mean_NResp.MI_SEM=splitapply(sem,Data_NResp.MI,G);
mean_NResp.Amp=splitapply(omean,Data_NResp.PeaksMax_mean,G);
mean_NResp.Amp_SEM=splitapply(sem,Data_NResp.PeaksMax_mean,G);
mean_NResp.Latency=splitapply(omean,Data_NResp.PeaksLatency_mean,G);
mean_NResp.Latency_SEM=splitapply(sem,Data_NResp.PeaksLatency_mean,G);
mean_NResp.PrctStimActive=splitapply(omean,Data_NResp.PrctStimActive,G);
mean_NResp.PrctStimActive_SEM=splitapply(sem,Data_NResp.PrctStimActive,G);
mean_NResp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp.allstimsROIs_mean),G);
mean_NResp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp.allstimsROIs_mean),G);

% Determine an MI threshold for responsivity
% Plot cells identified as responsive by bootstrapping stratified by MI
dataset = Data(Data.Timepoint == "BL" & Data.ResponsivityStim == 1,:);
dataset_1=dataset(dataset.MI<0.5,:); %There are 6/558 responders w/MI<0.5
dataset_2=dataset(dataset.MI>=0.5 & dataset.MI<1,:); %There are 18/558 responders w/MI 0.5-1
dataset_3=dataset(dataset.MI>=1 & dataset.MI<2,:); %There are 91/558 responders w/MI 1-2

%Plot responsive cells MI <0.5
figure_meanEvokedResp1 = figure;
axes_meanEvokedResp1 = axes('Parent',figure_meanEvokedResp1);
hold(axes_meanEvokedResp1,'on');

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

mean_trace=mean(cell2mat(dataset_1.allstimsROIs_mean),1);
greybar=NaN(size(mean_trace,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(dataset_1,1)
    plot(x, dataset_1.allstimsROIs_mean{i})
    hold on
%     pause
end

plot(x, mean_trace,'-k','LineWidth',4)
hold off

% Add the x-axis labels and title
set(axes_meanEvokedResp1, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of MI<0.5 responsive ROIs at baseline'); % Choice of the name of the plot HERE
ylim([-1 14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

%Plot for MI 0.5 to 1
figure_meanEvokedResp2 = figure;
axes_meanEvokedResp2 = axes('Parent',figure_meanEvokedResp2);
hold(axes_meanEvokedResp2,'on');

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

mean_trace=mean(cell2mat(dataset_2.allstimsROIs_mean),1);
greybar=NaN(size(mean_trace,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(dataset_2,1)
    plot(x, dataset_2.allstimsROIs_mean{i})
    hold on
%     pause
end

plot(x, mean_trace,'-k','LineWidth',4)
hold off

% Add the x-axis labels and title
set(axes_meanEvokedResp2, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of MI 0.5-1 responsive ROIs at baseline'); % Choice of the name of the plot HERE
ylim([-1 14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

%Plot for MI 1 to 2
figure_meanEvokedResp3 = figure;
axes_meanEvokedResp3 = axes('Parent',figure_meanEvokedResp3);
hold(axes_meanEvokedResp3,'on');

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

mean_trace=mean(cell2mat(dataset_3.allstimsROIs_mean),1);
greybar=NaN(size(mean_trace,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(dataset_3,1)
    plot(x, dataset_3.allstimsROIs_mean{i})
    hold on
%     pause
end

plot(x, mean_trace,'-k','LineWidth',4)
hold off

% Add the x-axis labels and title
set(axes_meanEvokedResp3, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of MI 1-2 responsive ROIs at baseline'); % Choice of the name of the plot HERE
ylim([-1 14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

clear xLabels mean_trace dataset dataset_1 dataset_2 dataset_3 greybar i x
clearvars -except Data Data_ByCell Data_ByCell_NoTrim Data_ByCell_wCorr...
    data_core Data_NoTrim Spont_Data_ByCell Spont_Data_ByCell_wCorr...
    mean_ByMovie Spont_mean_ByMovie
%% Section 2 (Required): Add MI to criteria for determining responsive cells
% Responsivity to whisker stimulation was determined as significant
% correlation between the calcium fluoresence trace
% ('signal_percentile_stim<0.01) and MI>=1. 

% Use a MI of >1. If you want to use a different MI, change the value on
% line 294, or remove '&& Data.MI(i)<1' to analyze the data without any MI
% criteria.
for i=1:size(Data,1)
    if Data.ResponsivityStim(i) == 1 && Data.MI(i)<1 %If a responsive cell has MI<1, set to non-responder
        Data.ResponsivityStim(i) = 0;
    else
    end
end

%Separate cells out by responders and non-responders in each barrel imaged
%at baseline
C1_BL = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL",:);
C1_BL_resp = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 1,:);
C1_BL_non_resp = Data(Data.Barrel_imaged == "C1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 0,:);
D1_BL = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL",:);

%Get means for each barrel/stim/timepoint from individual mice using shuffled stimulus responsivity criteria
omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);

%First calculate the percentage of whisker responsive cells within a mouse
T1=Data(:,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim','Gender'});
[G,mean_ByMouse]=findgroups(T1);
C=groupcounts(T1,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim','Gender'});

mean_ByMouse.NumRespCells=splitapply(osum,Data.ResponsivityStim,G);
mean_ByMouse.TotalCells=C.GroupCount;
mean_ByMouse.PercentResponsive=(mean_ByMouse.NumRespCells./mean_ByMouse.TotalCells)*100;

mean_ByMouse.MI_NumRespCells=splitapply(osum,Data.MI_ResponsivityStim,G);
mean_ByMouse.MI_PercentResponsive=(mean_ByMouse.MI_NumRespCells./mean_ByMouse.TotalCells)*100;

%Then summarize by group to get mean with n=# of mice
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_ByGroup]=findgroups(T1);

mean_ByGroup.PercentResponsive=splitapply(omean,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.PercentResponsive_SEM=splitapply(sem,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.MI_PercentResponsive=splitapply(omean,mean_ByMouse.MI_PercentResponsive,G);
mean_ByGroup.MI_PercentResponsive_SEM=splitapply(sem,mean_ByMouse.MI_PercentResponsive,G);

% Pool all cells by group for responsive cells and non-responsive cells
% First the responsive cells
Data_Resp=Data(Data.ResponsivityStim==1,:);
T1=Data_Resp(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_Resp]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,1,'omitnan')/sqrt(size(x,1)));

mean_Resp.TotalCells=C.GroupCount;
mean_Resp.AUC=splitapply(omean,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.AUC_SEM=splitapply(sem,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.MI=splitapply(omean,Data_Resp.MI,G);
mean_Resp.MI_SEM=splitapply(sem,Data_Resp.MI,G);
mean_Resp.Amp=splitapply(omean,Data_Resp.PeaksMax_mean,G);
mean_Resp.Amp_SEM=splitapply(sem,Data_Resp.PeaksMax_mean,G);
mean_Resp.Latency=splitapply(omean,Data_Resp.PeaksLatency_mean,G);
mean_Resp.Latency_SEM=splitapply(sem,Data_Resp.PeaksLatency_mean,G);
mean_Resp.PrctStimActive=splitapply(omean,Data_Resp.PrctStimActive,G);
mean_Resp.PrctStimActive_SEM=splitapply(sem,Data_Resp.PrctStimActive,G);
mean_Resp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp.allstimsROIs_mean),G);
mean_Resp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp.allstimsROIs_mean),G);

%Then the non-responsive cells
Data_NResp=Data(Data.ResponsivityStim==0,:);
T1=Data_NResp(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_NResp]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});

mean_NResp.TotalCells=C.GroupCount;
mean_NResp.AUC=splitapply(omean,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.AUC_SEM=splitapply(sem,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.MI=splitapply(omean,Data_NResp.MI,G);
mean_NResp.MI_SEM=splitapply(sem,Data_NResp.MI,G);
mean_NResp.Amp=splitapply(omean,Data_NResp.PeaksMax_mean,G);
mean_NResp.Amp_SEM=splitapply(sem,Data_NResp.PeaksMax_mean,G);
mean_NResp.Latency=splitapply(omean,Data_NResp.PeaksLatency_mean,G);
mean_NResp.Latency_SEM=splitapply(sem,Data_NResp.PeaksLatency_mean,G);
mean_NResp.PrctStimActive=splitapply(omean,Data_NResp.PrctStimActive,G);
mean_NResp.PrctStimActive_SEM=splitapply(sem,Data_NResp.PrctStimActive,G);
mean_NResp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp.allstimsROIs_mean),G);
mean_NResp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp.allstimsROIs_mean),G);

% For spontaneous data, pool all cells by group
T1=Spont_Data_ByCell(:,{'Timepoint','Barrel_imaged'});
[G,spont_mean_allcells]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged'});

spont_mean_allcells.TotalCells=C.GroupCount;
spont_mean_allcells.AUC=splitapply(omean,Spont_Data_ByCell.AUC,G);
spont_mean_allcells.AUC_SEM=splitapply(sem,Spont_Data_ByCell.AUC,G);
spont_mean_allcells.mean_Z=splitapply(omean,Spont_Data_ByCell.mean_Z,G);
spont_mean_allcells.mean_Z_SEM=splitapply(sem,Spont_Data_ByCell.mean_Z,G);
spont_mean_allcells.PeaksAmp_mean=splitapply(omean,Spont_Data_ByCell.PeaksAmp_mean,G);
spont_mean_allcells.PeaksAmp_mean_SEM=splitapply(sem,Spont_Data_ByCell.PeaksAmp_mean,G);
spont_mean_allcells.PeakFreq=splitapply(omean,Spont_Data_ByCell.PeakFreq,G);
spont_mean_allcells.PeakFreq_SEM=splitapply(sem,Spont_Data_ByCell.PeakFreq,G);

% For No Trim Control data, Get means for each barrel/stim/timepoint from individual mice
T1=Data_NoTrim(:,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_ByMouse_NoTrim]=findgroups(T1);
C=groupcounts(T1,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim'});

mean_ByMouse_NoTrim.NumRespCells=splitapply(osum,Data_NoTrim.Responsivity,G);
mean_ByMouse_NoTrim.TotalCells=C.GroupCount;
mean_ByMouse_NoTrim.PercentResponsive=(mean_ByMouse_NoTrim.NumRespCells./mean_ByMouse_NoTrim.TotalCells)*100;

%Then summarize by group to get mean with n=# of mice
T1=mean_ByMouse_NoTrim(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_ByGroup_NoTrim]=findgroups(T1);

mean_ByGroup_NoTrim.PercentResponsive=splitapply(omean,mean_ByMouse_NoTrim.PercentResponsive,G);
mean_ByGroup_NoTrim.PercentResponsive_SEM=splitapply(sem,mean_ByMouse_NoTrim.PercentResponsive,G);

% Pool all cells by group for responsive cells and non-responsive cells
Data_Resp_NoTrim=Data_NoTrim(Data_NoTrim.Responsivity==1,:);
T1=Data_Resp_NoTrim(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_Resp_NoTrim]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});

mean_Resp_NoTrim.TotalCells=C.GroupCount;
mean_Resp_NoTrim.AUC=splitapply(omean,Data_Resp_NoTrim.allstimsallROIs_AUC,G);
mean_Resp_NoTrim.AUC_SEM=splitapply(sem,Data_Resp_NoTrim.allstimsallROIs_AUC,G);
mean_Resp_NoTrim.MI=splitapply(omean,Data_Resp_NoTrim.MI,G);
mean_Resp_NoTrim.MI_SEM=splitapply(sem,Data_Resp_NoTrim.MI,G);
mean_Resp_NoTrim.Amp=splitapply(omean,Data_Resp_NoTrim.PeaksMax_mean,G);
mean_Resp_NoTrim.Amp_SEM=splitapply(sem,Data_Resp_NoTrim.PeaksMax_mean,G);
mean_Resp_NoTrim.Latency=splitapply(omean,Data_Resp_NoTrim.PeaksLatency_mean,G);
mean_Resp_NoTrim.Latency_SEM=splitapply(sem,Data_Resp_NoTrim.PeaksLatency_mean,G);
mean_Resp_NoTrim.PrctStimActive=splitapply(omean,Data_Resp_NoTrim.PrctStimActive,G);
mean_Resp_NoTrim.PrctStimActive_SEM=splitapply(sem,Data_Resp_NoTrim.PrctStimActive,G);
mean_Resp_NoTrim.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp_NoTrim.allstimsROIs_mean),G);
mean_Resp_NoTrim.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp_NoTrim.allstimsROIs_mean),G);

Data_NResp_NoTrim=Data_NoTrim(Data_NoTrim.Responsivity==0,:);
T1=Data_NResp_NoTrim(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G,mean_NResp_NoTrim]=findgroups(T1);
C=groupcounts(T1,{'Timepoint','Barrel_imaged','Whisker_stim'});

mean_NResp_NoTrim.TotalCells=C.GroupCount;
mean_NResp_NoTrim.AUC=splitapply(omean,Data_NResp_NoTrim.allstimsallROIs_AUC,G);
mean_NResp_NoTrim.AUC_SEM=splitapply(sem,Data_NResp_NoTrim.allstimsallROIs_AUC,G);
mean_NResp_NoTrim.MI=splitapply(omean,Data_NResp_NoTrim.MI,G);
mean_NResp_NoTrim.MI_SEM=splitapply(sem,Data_NResp_NoTrim.MI,G);
mean_NResp_NoTrim.Amp=splitapply(omean,Data_NResp_NoTrim.PeaksMax_mean,G);
mean_NResp_NoTrim.Amp_SEM=splitapply(sem,Data_NResp_NoTrim.PeaksMax_mean,G);
mean_NResp_NoTrim.Latency=splitapply(omean,Data_NResp_NoTrim.PeaksLatency_mean,G);
mean_NResp_NoTrim.Latency_SEM=splitapply(sem,Data_NResp_NoTrim.PeaksLatency_mean,G);
mean_NResp_NoTrim.PrctStimActive=splitapply(omean,Data_NResp_NoTrim.PrctStimActive,G);
mean_NResp_NoTrim.PrctStimActive_SEM=splitapply(sem,Data_NResp_NoTrim.PrctStimActive,G);
mean_NResp_NoTrim.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp_NoTrim.allstimsROIs_mean),G);
mean_NResp_NoTrim.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp_NoTrim.allstimsROIs_mean),G);

clear i T1 G C omean omean2 osum sem sem2

%% Section 3 (Required): Assign whisker responses, calculate SI, and switch categories
    Data = sortrows(Data, "Index");  

    % Selection of barrel, timepoint & of mice whose the 4 whiskers have been stimulated
    C1_BL_D46 = Data(Data.Barrel_imaged == "C1" & (Data.Timepoint == "BL"|Data.Timepoint =="Day46" ) & (Data.MouseID == "m582"|Data.MouseID == "m584"|Data.MouseID == "m585"|Data.MouseID == "m604"|Data.MouseID == "m605"|Data.MouseID == "m607"),:);

    % Creation of the dataset whose each row represent a unique PV cell
    % (this is why it is important that the order is BL-D46 and B2-C1-D1-E3)
    % --> manually check the order of the rows !!!! 

    dataset = C1_BL_D46;

    New_data = table;
    for row=1:8:height(dataset)
        New_data.MouseID(row) = dataset.MouseID(row);
        New_data.Group(row) = dataset.Group(row);
        New_data.Barrel_imaged(row) = dataset.Barrel_imaged(row);
        New_data.Index(row) = dataset.Index(row);

        New_data.AUC_BL_B2(row) = dataset.allstimsallROIs_AUC(row);
        New_data.AUC_BL_C1(row) = dataset.allstimsallROIs_AUC(row+1);
        New_data.AUC_BL_D1(row) = dataset.allstimsallROIs_AUC(row+2);
        New_data.AUC_BL_E3(row) = dataset.allstimsallROIs_AUC(row+3);
        New_data.AUC_D46_B2(row) = dataset.allstimsallROIs_AUC(row+4);
        New_data.AUC_D46_C1(row) = dataset.allstimsallROIs_AUC(row+5);
        New_data.AUC_D46_D1(row) = dataset.allstimsallROIs_AUC(row+6);
        New_data.AUC_D46_E3(row) = dataset.allstimsallROIs_AUC(row+7);
        New_data.MI_BL_B2(row) = dataset.MI(row);
        New_data.MI_BL_C1(row) = dataset.MI(row+1);
        New_data.MI_BL_D1(row) = dataset.MI(row+2);
        New_data.MI_BL_E3(row) = dataset.MI(row+3);
        New_data.MI_D46_B2(row) = dataset.MI(row+4);
        New_data.MI_D46_C1(row) = dataset.MI(row+5);
        New_data.MI_D46_D1(row) = dataset.MI(row+6);
        New_data.MI_D46_E3(row) = dataset.MI(row+7);
        New_data.Responsiveness_BL_B2(row) = dataset.ResponsivityStim(row);
        New_data.Responsiveness_BL_C1(row) = dataset.ResponsivityStim(row+1);
        New_data.Responsiveness_BL_D1(row) = dataset.ResponsivityStim(row+2);
        New_data.Responsiveness_BL_E3(row) = dataset.ResponsivityStim(row+3);
        New_data.Responsiveness_D46_B2(row) = dataset.ResponsivityStim(row+4);
        New_data.Responsiveness_D46_C1(row) = dataset.ResponsivityStim(row+5);
        New_data.Responsiveness_D46_D1(row) = dataset.ResponsivityStim(row+6);
        New_data.Responsiveness_D46_E3(row) = dataset.ResponsivityStim(row+7);

    end

    clear row  dataset

    % Suppresion of the empty rows (we go through the list in oposite order because we use the index of the rows, and we suppress them)
    for row = height(New_data):-1:1
        if New_data.Index(row) == 0
            New_data(row,:)= [];
        end
    end
    clear row


    C1_BL_D46 = New_data; 
    clear New_data

% Now for the C1 barrel, only C1 and D1 whisker stimulations @ BL and D46
    C1_BL_D46_all  = Data(Data.Barrel_imaged == "C1" & (Data.Timepoint == "BL"|Data.Timepoint =="Day46" ) & (Data.Whisker_stim =='Stm_C1'|Data.Whisker_stim=='Stm_D1'),:);
    dataset = C1_BL_D46_all;

    New_data = table;
    for row=1:4:height(dataset)
        New_data.MouseID(row) = dataset.MouseID(row);
        New_data.Group(row) = dataset.Group(row);
        New_data.Barrel_imaged(row) = dataset.Barrel_imaged(row);
        New_data.Index(row) = dataset.Index(row);

        New_data.AUC_BL_C1(row) = dataset.allstimsallROIs_AUC(row);
        New_data.AUC_BL_D1(row) = dataset.allstimsallROIs_AUC(row+1);
        New_data.AUC_D46_C1(row) = dataset.allstimsallROIs_AUC(row+2);
        New_data.AUC_D46_D1(row) = dataset.allstimsallROIs_AUC(row+3);
        New_data.MI_BL_C1(row) = dataset.MI(row);
        New_data.MI_BL_D1(row) = dataset.MI(row+1);
        New_data.MI_D46_C1(row) = dataset.MI(row+2);
        New_data.MI_D46_D1(row) = dataset.MI(row+3);
        New_data.Responsiveness_BL_C1(row) = dataset.ResponsivityStim(row);
        New_data.Responsiveness_BL_D1(row) = dataset.ResponsivityStim(row+1);
        New_data.Responsiveness_D46_C1(row) = dataset.ResponsivityStim(row+2);
        New_data.Responsiveness_D46_D1(row) = dataset.ResponsivityStim(row+3);  
    end

    clear row  dataset

    % Suppresion of the empty rows (we go through the list in oposite order because we use the index of the rows, and we suppress them)
    for row = height(New_data):-1:1
        if New_data.Index(row) == 0
            New_data(row,:)= [];
        end
    end
    clear row


    C1_BL_D46_all = New_data; 
    clear New_data

% Now for the D1 barrel, C1 and D1 whisker stimulations @ BL and D46
    D1_BL_D46  = Data(Data.Barrel_imaged == "D1" & (Data.Timepoint == "BL"|Data.Timepoint =="Day46" ),:);
    dataset = D1_BL_D46;

    New_data = table;
    for row=1:4:height(dataset)
        New_data.MouseID(row) = dataset.MouseID(row);
        New_data.Group(row) = dataset.Group(row);
        New_data.Barrel_imaged(row) = dataset.Barrel_imaged(row);
        New_data.Index(row) = dataset.Index(row);

        New_data.AUC_BL_C1(row) = dataset.allstimsallROIs_AUC(row);
        New_data.AUC_BL_D1(row) = dataset.allstimsallROIs_AUC(row+1);
        New_data.AUC_D46_C1(row) = dataset.allstimsallROIs_AUC(row+2);
        New_data.AUC_D46_D1(row) = dataset.allstimsallROIs_AUC(row+3);
        New_data.MI_BL_C1(row) = dataset.MI(row);
        New_data.MI_BL_D1(row) = dataset.MI(row+1);
        New_data.MI_D46_C1(row) = dataset.MI(row+2);
        New_data.MI_D46_D1(row) = dataset.MI(row+3);
        New_data.Responsiveness_BL_C1(row) = dataset.ResponsivityStim(row);
        New_data.Responsiveness_BL_D1(row) = dataset.ResponsivityStim(row+1);
        New_data.Responsiveness_D46_C1(row) = dataset.ResponsivityStim(row+2);
        New_data.Responsiveness_D46_D1(row) = dataset.ResponsivityStim(row+3);  
    end

    clear row  dataset

    % Suppresion of the empty rows (we go through the list in oposite order because we use the index of the rows, and we suppress them)
    for row = height(New_data):-1:1
        if New_data.Index(row) == 0
            New_data(row,:)= [];
        end
    end
    clear row


    D1_BL_D46 = New_data; 
    clear New_data
    
% C1 barrel, changes in responsivity to D1 stim over time
% (this is why it is important that the order is BL-Day3-Day6-Day12-Day18-Day46)
% --> manually check the order of the rows !!!!
C1_Stm_D1 = Data(Data.Barrel_imaged == "C1" & Data.Whisker_stim == "Stm_D1",:);

dataset = C1_Stm_D1;

New_data = table;
for row=1:6:height(dataset)
    New_data.MouseID(row) = dataset.MouseID(row);
    New_data.Group(row) = dataset.Group(row);
    New_data.Barrel_imaged(row) = dataset.Barrel_imaged(row);
    New_data.Whisker_stim(row) = dataset.Whisker_stim(row);
    New_data.Index(row) = dataset.Index(row);


    New_data.MI_BL_D1(row) = dataset.MI(row);
    New_data.MI_D3_D1(row) = dataset.MI(row+1);
    New_data.MI_D6_D1(row) = dataset.MI(row+2);
    New_data.MI_D12_D1(row) = dataset.MI(row+3);
    New_data.MI_D18_D1(row) = dataset.MI(row+4);
    New_data.MI_D46_D1(row) = dataset.MI(row+5);
    New_data.Responsiveness_BL_D1(row) = dataset.ResponsivityStim(row);
    New_data.Responsiveness_D3_D1(row) = dataset.ResponsivityStim(row+1);
    New_data.Responsiveness_D6_D1(row) = dataset.ResponsivityStim(row+2);
    New_data.Responsiveness_D12_D1(row) = dataset.ResponsivityStim(row+3);
    New_data.Responsiveness_D18_D1(row) = dataset.ResponsivityStim(row+4);
    New_data.Responsiveness_D46_D1(row) = dataset.ResponsivityStim(row+5);
end

clear row  dataset

% Suppresion of the empty rows (we go through the list in oposite order because we use the index of the rows, and we suppress them)
for row = height(New_data):-1:1
    if New_data.Index(row) == 0
        New_data(row,:)= [];
    end
end
clear row

C1_Stm_D1 = New_data; 
clear New_data

% D1 barrel, changes in responsivity to D1 stim over time
% (this is why it is important that the order is BL-Day3-Day6-Day12-Day18-Day46)
% --> manually check the order of the rows !!!!
D1_Stm_D1 = Data(Data.Barrel_imaged == "D1" & Data.Whisker_stim == "Stm_D1",:);

dataset = D1_Stm_D1;

New_data = table;
for row=1:6:height(dataset)
    New_data.MouseID(row) = dataset.MouseID(row);
    New_data.Group(row) = dataset.Group(row);
    New_data.Barrel_imaged(row) = dataset.Barrel_imaged(row);
    New_data.Whisker_stim(row) = dataset.Whisker_stim(row);
    New_data.Index(row) = dataset.Index(row);


    New_data.MI_BL_D1(row) = dataset.MI(row);
    New_data.MI_D3_D1(row) = dataset.MI(row+1);
    New_data.MI_D6_D1(row) = dataset.MI(row+2);
    New_data.MI_D12_D1(row) = dataset.MI(row+3);
    New_data.MI_D18_D1(row) = dataset.MI(row+4);
    New_data.MI_D46_D1(row) = dataset.MI(row+5);
    New_data.Responsiveness_BL_D1(row) = dataset.ResponsivityStim(row);
    New_data.Responsiveness_D3_D1(row) = dataset.ResponsivityStim(row+1);
    New_data.Responsiveness_D6_D1(row) = dataset.ResponsivityStim(row+2);
    New_data.Responsiveness_D12_D1(row) = dataset.ResponsivityStim(row+3);
    New_data.Responsiveness_D18_D1(row) = dataset.ResponsivityStim(row+4);
    New_data.Responsiveness_D46_D1(row) = dataset.ResponsivityStim(row+5);
end

clear row  dataset

% Suppresion of the empty rows (we go through the list in oposite order because we use the index of the rows, and we suppress them)
for row = height(New_data):-1:1
    if New_data.Index(row) == 0
        New_data(row,:)= [];
    end
end
clear row

D1_Stm_D1 = New_data; 
clear New_data


% --------------------Calculate Selectivity Index---------------------
    
% Selectivity Index in C1_BL_D46_all at BL and at D46 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = C1_BL_D46_all;
% Set negative MI & AUC to 0 for calculation of SI or TCW
for row=1:height(dataset)
    for col=5:12
        if dataset{row,col} < 0
            dataset{row,col}= 0;
        end
    end
end
clear row col

for row=1:height(dataset)
    % Computation of SI: (AUC_C1 - AUC_D1) / (AUC_C1 + AUC_D1)
    % SI belongs to [-1,1]
    % SI = 0 means no preference for either C1 or D1 (same response)
    % SI = -1 means total selectivity for D1
    % SI = 1 means total selectivity for C1
    dataset.SI_BL(row) = (dataset.AUC_BL_C1(row) - dataset.AUC_BL_D1(row)) / (dataset.AUC_BL_C1(row) + dataset.AUC_BL_D1(row));
    dataset.SI_D46(row) = (dataset.AUC_D46_C1(row) - dataset.AUC_D46_D1(row)) / (dataset.AUC_D46_D1(row) + dataset.AUC_D46_C1(row));
    dataset.SI_MI_BL(row) = (dataset.MI_BL_C1(row) - dataset.MI_BL_D1(row)) / (dataset.MI_BL_C1(row) + dataset.MI_BL_D1(row));
    dataset.SI_MI_D46(row) = (dataset.MI_D46_C1(row) - dataset.MI_D46_D1(row)) / (dataset.MI_D46_D1(row) + dataset.MI_D46_C1(row));
end

%Set SI=NaN (meaning no response to C1 or D1 whisker) to 0
dataset.SI_BL(isnan(dataset.SI_BL))=0;
dataset.SI_D46(isnan(dataset.SI_D46))=0;
dataset.SI_MI_BL(isnan(dataset.SI_MI_BL))=0;
dataset.SI_MI_D46(isnan(dataset.SI_MI_D46))=0;

C1_BL_D46_all = dataset;
clear dataset

% Selectivity Index in D1_BL_D46 at BL and at D46 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset = D1_BL_D46;

% Set negative MI & AUC to 0 for calculation of SI or TCW
for row=1:height(dataset)
    for col=5:12
        if dataset{row,col} < 0
            dataset{row,col}= 0;
        end
    end
end
clear row col

for row=1:height(dataset)
    % Computation of SI: (AUC_C1 - AUC_D1) / (AUC_C1 + AUC_D1)
    % SI belongs to [-1,1]
    % SI = 0 means no preference for either C1 or D1 (same response)
    % SI = -1 means total selectivity for D1
    % SI = 1 means total selectivity for C1
    dataset.SI_BL(row) = (dataset.AUC_BL_C1(row) - dataset.AUC_BL_D1(row)) / (dataset.AUC_BL_C1(row) + dataset.AUC_BL_D1(row));
    dataset.SI_D46(row) = (dataset.AUC_D46_C1(row) - dataset.AUC_D46_D1(row)) / (dataset.AUC_D46_D1(row) + dataset.AUC_D46_C1(row));
    dataset.SI_MI_BL(row) = (dataset.MI_BL_C1(row) - dataset.MI_BL_D1(row)) / (dataset.MI_BL_C1(row) + dataset.MI_BL_D1(row));
    dataset.SI_MI_D46(row) = (dataset.MI_D46_C1(row) - dataset.MI_D46_D1(row)) / (dataset.MI_D46_D1(row) + dataset.MI_D46_C1(row));
end

%Set SI=NaN (meaning no response to C1 or D1 whisker) to 0
dataset.SI_BL(isnan(dataset.SI_BL))=0;
dataset.SI_D46(isnan(dataset.SI_D46))=0;
dataset.SI_MI_BL(isnan(dataset.SI_MI_BL))=0;
dataset.SI_MI_D46(isnan(dataset.SI_MI_D46))=0;

D1_BL_D46 = dataset;
clear dataset row i


% -----------Assign switch categories BL to D12, BL to D46---------------
%C1 barrel with D1 stimulation
dataset = C1_Stm_D1;
%------ BL to D12 --------
for row=1:height(dataset)
    if double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 1
        dataset.Switch_BL_to_D12(row) = {'NR to R'}; 
    elseif double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == -1
        dataset.Switch_BL_to_D12(row) = {'R to NR'};
    elseif double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 0
        if double(string(dataset.Responsiveness_BL_D1(row))) == 0
            dataset.Switch_BL_to_D12(row) = {'NR to NR'};
        elseif double(string(dataset.Responsiveness_BL_D1(row))) == 1
            dataset.Switch_BL_to_D12(row) = {'R to R'};
        end
    end
end
clear row
dataset.Switch_BL_to_D12=categorical(dataset.Switch_BL_to_D12);
% Category order: NR to NR, NR to R, R to NR, R to R
%------ BL to D46 --------
for row=1:height(dataset)
    if double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 1
        dataset.Switch_BL_to_D46(row) = {'NR to R'}; 
    elseif double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == -1
        dataset.Switch_BL_to_D46(row) = {'R to NR'};
    elseif double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 0
        if double(string(dataset.Responsiveness_BL_D1(row))) == 0
            dataset.Switch_BL_to_D46(row) = {'NR to NR'};
        elseif double(string(dataset.Responsiveness_BL_D1(row))) == 1
            dataset.Switch_BL_to_D46(row) = {'R to R'};
        end
    end
end
dataset.Switch_BL_to_D46=categorical(dataset.Switch_BL_to_D46);
% Category order: NR to NR, NR to R, R to NR, R to R

% Add in SI @ baseline
dataset.SI_MI_BL=C1_BL_D46_all.SI_MI_BL;

C1_Stm_D1 = dataset;

%D1 barrel with D1 stimulation
dataset = D1_Stm_D1;
%------ BL to D12 --------
for row=1:height(dataset)
    if double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 1
        dataset.Switch_BL_to_D12(row) = {'NR to R'}; 
    elseif double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == -1
        dataset.Switch_BL_to_D12(row) = {'R to NR'};
    elseif double(string(dataset.Responsiveness_D12_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 0
        if double(string(dataset.Responsiveness_BL_D1(row))) == 0
            dataset.Switch_BL_to_D12(row) = {'NR to NR'};
        elseif double(string(dataset.Responsiveness_BL_D1(row))) == 1
            dataset.Switch_BL_to_D12(row) = {'R to R'};
        end
    end
end
clear row
dataset.Switch_BL_to_D12=categorical(dataset.Switch_BL_to_D12);
% Category order: NR to NR, NR to R, R to NR, R to R

%------ BL to D46 --------
for row=1:height(dataset)
    if double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 1
        dataset.Switch_BL_to_D46(row) = {'NR to R'}; 
    elseif double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == -1
        dataset.Switch_BL_to_D46(row) = {'R to NR'};
    elseif double(string(dataset.Responsiveness_D46_D1(row))) - double(string(dataset.Responsiveness_BL_D1(row))) == 0
        if double(string(dataset.Responsiveness_BL_D1(row))) == 0
            dataset.Switch_BL_to_D46(row) = {'NR to NR'};
        elseif double(string(dataset.Responsiveness_BL_D1(row))) == 1
            dataset.Switch_BL_to_D46(row) = {'R to R'};
        end
    end
end
dataset.Switch_BL_to_D46=categorical(dataset.Switch_BL_to_D46);
% Category order: NR to NR, NR to R, R to NR, R to R

% Add in SI @ baseline
dataset.SI_MI_BL=D1_BL_D46.SI_MI_BL;

D1_Stm_D1 = dataset;

%Add switch categories to C1 barrel C1 stimulation:
C1_Stm_C1_BL= Data(Data.Barrel_imaged == "C1" & Data.Whisker_stim == "Stm_C1" & Data.Timepoint=='BL',:);
C1_Stm_C1_BL.Switch_BL_to_D12=C1_Stm_D1.Switch_BL_to_D12;
C1_Stm_C1_BL.Switch_BL_to_D46=C1_Stm_D1.Switch_BL_to_D46;
C1_Stm_C1_BL.SI_MI_BL=C1_Stm_D1.SI_MI_BL;

%Add switch categories to D1 barrel C1 stimulation:
D1_Stm_C1_BL= Data(Data.Barrel_imaged == "D1" & Data.Whisker_stim == "Stm_C1" & Data.Timepoint=='BL',:);
D1_Stm_C1_BL.Switch_BL_to_D12=D1_Stm_D1.Switch_BL_to_D12;
D1_Stm_C1_BL.Switch_BL_to_D46=D1_Stm_D1.Switch_BL_to_D46;
D1_Stm_C1_BL.SI_MI_BL=D1_Stm_D1.SI_MI_BL;

clear dataset row
%% Fig. 1d: Plot the % responsive cells to each whisker in the C1 barrel at baseline
dataset = mean_ByGroup(mean_ByGroup.Timepoint=='BL' & mean_ByGroup.Barrel_imaged=='C1',:); 

mean_C1 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_C1'); std_C1 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_C1');
mean_D1 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_D1'); std_D1 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_D1');
mean_E3 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_E3'); std_E3 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_E3');
mean_B2 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_B2'); std_B2 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_B2');

means = [mean_C1; mean_D1; mean_B2; mean_E3;]; stdErrors = [std_C1; std_D1; std_B2; std_E3;];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1', 'B2', 'E3'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_PrctResp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none'); % parameters to modify for the design (if needed)
dataset = mean_ByMouse(mean_ByMouse.Timepoint=='BL' & mean_ByMouse.Barrel_imaged=='C1',:);
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:size(whiskers,1)
    scatter(k*(ones(length(dataset.PercentResponsive(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.PercentResponsive(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('% Whisker Responsive Cells in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability

%Get statistics
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Whisker_stim + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells,'FitMethod','REMPL')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means
clear adj_ci_cvrg adj_p h crit_p p tbl stats c colors glme k normality whiskers dataset
%% Fig. 1e: Plot the mean-evoked traces from all cells to each whisker in the C1 barrel @ baseline
dataset = C1_BL; 


Stm_C1 = []; Stm_D1 = []; Stm_B2 = []; Stm_E3 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_B2" 
        Stm_B2 = [Stm_B2 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_E3" 
        Stm_E3 = [Stm_E3 dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(cell2mat(Stm_C1'), "omitnan"); std_C1 = std(cell2mat(Stm_C1'), "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(cell2mat(Stm_D1'), "omitnan"); std_D1 = std(cell2mat(Stm_D1'), "omitnan")/sqrt(width(Stm_D1));
mean_E3 = mean(cell2mat(Stm_E3'), "omitnan"); std_E3 = std(cell2mat(Stm_E3'), "omitnan")/sqrt(width(Stm_E3));
mean_B2 = mean(cell2mat(Stm_B2'), "omitnan"); std_B2 = std(cell2mat(Stm_B2'), "omitnan")/sqrt(width(Stm_B2));
clear Stm_C1 Stm_D1 Stm_B2 Stm_E3

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_C1', mean_D1', mean_B2', mean_E3']; stdErrors = [std_C1', std_D1', std_B2', std_E3'];

figure_meanEvokedAll = figure;
axes_meanEvokedAll = axes('Parent',figure_meanEvokedAll);
hold(axes_meanEvokedAll,'on');
style={'-m','-b','-g','-y'}; %Define colors for each whisker stimulated

%Plot the stimulus epoch as a grey area
greybar=NaN(size(mean_C1,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');
hold on

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of all ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])


clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means dataset greybar i style x
 
%% Fig. 1f: Plot the mean-evoked traces from responders to each whisker in the C1 barrel
dataset = C1_BL_resp; 

Stm_C1 = []; Stm_D1 = []; Stm_B2 = []; Stm_E3 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_B2" 
        Stm_B2 = [Stm_B2 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_E3" 
        Stm_E3 = [Stm_E3 dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(cell2mat(Stm_C1'), "omitnan"); std_C1 = std(cell2mat(Stm_C1'), "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(cell2mat(Stm_D1'), "omitnan"); std_D1 = std(cell2mat(Stm_D1'), "omitnan")/sqrt(width(Stm_D1));
mean_E3 = mean(cell2mat(Stm_E3'), "omitnan"); std_E3 = std(cell2mat(Stm_E3'), "omitnan")/sqrt(width(Stm_E3));
mean_B2 = mean(cell2mat(Stm_B2'), "omitnan"); std_B2 = std(cell2mat(Stm_B2'), "omitnan")/sqrt(width(Stm_B2));
clear Stm_C1 Stm_D1 Stm_B2 Stm_E3

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_C1', mean_D1', mean_B2', mean_E3']; stdErrors = [std_C1', std_D1', std_B2', std_E3'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-m','-b','-g','-y'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_C1,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of responsive ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])


clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means dataset greybar i style x
%% Fig. 1g: Plot the mean-evoked traces from non-responders to each whisker in the C1 barrel at baseline
dataset = C1_BL_non_resp; 

Stm_C1 = []; Stm_D1 = []; Stm_B2 = []; Stm_E3 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_B2" 
        Stm_B2 = [Stm_B2 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_E3" 
        Stm_E3 = [Stm_E3 dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(cell2mat(Stm_C1'), "omitnan"); std_C1 = std(cell2mat(Stm_C1'), "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(cell2mat(Stm_D1'), "omitnan"); std_D1 = std(cell2mat(Stm_D1'), "omitnan")/sqrt(width(Stm_D1));
mean_E3 = mean(cell2mat(Stm_E3'), "omitnan"); std_E3 = std(cell2mat(Stm_E3'), "omitnan")/sqrt(width(Stm_E3));
mean_B2 = mean(cell2mat(Stm_B2'), "omitnan"); std_B2 = std(cell2mat(Stm_B2'), "omitnan")/sqrt(width(Stm_B2));
clear Stm_C1 Stm_D1 Stm_B2 Stm_E3

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_C1', mean_D1', mean_B2', mean_E3']; stdErrors = [std_C1', std_D1', std_B2', std_E3'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedNResp = figure;
axes_meanEvokedNResp = axes('Parent',figure_meanEvokedNResp);
hold(axes_meanEvokedNResp,'on');
style={'-m','-b','-g','-y'}; %Define colors for each whisker stimulated

%Plot the stimulus epoch as a grey area
greybar=NaN(size(mean_C1,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');
hold on

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedNResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of non-responsive ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means dataset greybar i style x
%% Fig. 1h: Plot the MI from Responders in the C1 barrel to all whiskers at baseline

% Choice of the dataset HERE
dataset = C1_BL_resp;

Stm_C1 = []; Stm_D1 = []; Stm_B2 = []; Stm_E3 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_B2" 
        Stm_B2 = [Stm_B2 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_E3" 
        Stm_E3 = [Stm_E3 dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(Stm_C1, "omitnan"); std_C1 = std(Stm_C1, "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(Stm_D1, "omitnan"); std_D1 = std(Stm_D1, "omitnan")/sqrt(width(Stm_D1));
mean_E3 = mean(Stm_E3, "omitnan"); std_E3 = std(Stm_E3, "omitnan")/sqrt(width(Stm_E3));
mean_B2 = mean(Stm_B2, "omitnan"); std_B2 = std(Stm_B2, "omitnan")/sqrt(width(Stm_B2));

xLabels = {'C1', 'D1', 'B2', 'E3'}; means = [mean_C1; mean_D1; mean_B2; mean_E3]; stdErrors = [std_C1; std_D1; std_B2; std_E3];
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];
% Plot the bar chart with error bars for the standard errors
figure_AUCresp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_AUCresp.CData(k,:) = colors(k,:);
end
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:size(whiskers,1)
    scatter(k*(ones(length(dataset.MI(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.MI(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Whisker stimulated'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of the responsive ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability

%Get statistics
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + (1 | MouseID)','FitMethod','REML') %random intercept
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear Stm_C1 Stm_D1 Stm_B2 Stm_E3
clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means
clear adj_ci_cvrg adj_p h crit_p colors lme k normality whiskers dataset
%% Fig. 1i: Plot the MI from Non-Responders in the C1 barrel to all whiskers at baseline

% Choice of the dataset HERE
dataset = C1_BL_non_resp;

Stm_C1 = []; Stm_D1 = []; Stm_B2 = []; Stm_E3 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_B2" 
        Stm_B2 = [Stm_B2 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_E3" 
        Stm_E3 = [Stm_E3 dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(Stm_C1, "omitnan"); std_C1 = std(Stm_C1, "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(Stm_D1, "omitnan"); std_D1 = std(Stm_D1, "omitnan")/sqrt(width(Stm_D1));
mean_E3 = mean(Stm_E3, "omitnan"); std_E3 = std(Stm_E3, "omitnan")/sqrt(width(Stm_E3));
mean_B2 = mean(Stm_B2, "omitnan"); std_B2 = std(Stm_B2, "omitnan")/sqrt(width(Stm_B2));
clear Stm_C1 Stm_D1 Stm_B2 Stm_E3

xLabels = {'C1', 'D1', 'B2', 'E3'}; means = [mean_C1; mean_D1; mean_B2; mean_E3]; stdErrors = [std_C1; std_D1; std_B2; std_E3];
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];
% Plot the bar chart with error bars for the standard errors
figure_AUCnresp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_AUCnresp.CData(k,:) = colors(k,:);
end
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:size(whiskers,1)
    scatter(k*(ones(length(dataset.MI(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.MI(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Whisker stimulated'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of the non-responsive ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability

%Get statistics
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means
clear adj_ci_cvrg adj_p h crit_p colors lme k normality whiskers dataset
%% Fig. 1j: Plot spatial distribution of responsive cells in C1 barrel at baseline
% Note data is mirrored across the x-axis to match the location of the cell in the imaging
% FOV, which has 0,0 at the NW corner, rather than MATLAB, which plots 0,0
% at the SW corner. Data are also flipped as the imaging FOV results the
% position as y,x rather than x,y.

%Plot C1 Whisker
dataset = C1_BL(C1_BL.Whisker_stim == "Stm_C1",:); 

f_C1_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color='m';
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('C1 Barrel - C1 Whisker Stimulation - Baseline');

%Plot D1 whisker
dataset = C1_BL(C1_BL.Whisker_stim == "Stm_D1",:); 

f_D1_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color='b';
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('C1 Barrel - D1 Whisker Stimulation - Baseline');

%Plot B2 Whisker
dataset = C1_BL(C1_BL.Whisker_stim == "Stm_B2",:); 

f_B2_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color='g';
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('C1 Barrel - B2 Whisker Stimulation - Baseline');

%Plot E3 whisker
dataset = C1_BL(C1_BL.Whisker_stim == "Stm_E3",:); 

f_E3_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color=[1 0.4 0.2];
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('C1 Barrel - E3 Whisker Stimulation - Baseline');

%Plot shapes for the legend
f_legend_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on
legend_rads=[1,5,20,40,80];
legend_centers=[10,256;50,256;100,256;200,256;400,256];

for i=1:size(legend_rads,2)
    drawcircle('Center',legend_centers(i,:),'Radius',legend_rads(i),'Color','k',...
    'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
    'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
drawrectangle('Position',[100,400,(50/0.825),10],'Color','k',...
    'FaceAlpha',1,'StripeColor','none','HandleVisibility','off',...
    'InteractionsAllowed','none','LineWidth',0.001);

title('spatial legend')

clear center color dataset i legend_centers legend_rads radius

%% Fig. S2a: Plot the % responsive cells to each whisker in the D1 barrel at baseline
dataset = mean_ByGroup(mean_ByGroup.Timepoint=='BL' & mean_ByGroup.Barrel_imaged=='D1',:); 

mean_C1 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_C1'); std_C1 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_C1');
mean_D1 = dataset.PercentResponsive(dataset.Whisker_stim=='Stm_D1'); std_D1 = dataset.PercentResponsive_SEM(dataset.Whisker_stim=='Stm_D1');

means = [mean_C1; mean_D1]; stdErrors = [std_C1; std_D1];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_PrctRespD1=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctRespD1.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none'); % parameters to modify for the design (if needed)
dataset = mean_ByMouse(mean_ByMouse.Timepoint=='BL' & mean_ByMouse.Barrel_imaged=='D1',:);
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:2 %Only C1 and D1 whiskers for D1 barrel
    scatter(k*(ones(length(dataset.PercentResponsive(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.PercentResponsive(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('% Whisker Responsive Cells in D1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
normality = lillietest(dataset.NumRespCells) %0, normal
dataset.Whisker_stim=removecats(dataset.Whisker_stim, {'Stm_B2','Stm_E3'});
glme = fitglme(dataset,'NumRespCells ~ Whisker_stim + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells,'FitMethod','REMPL')

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means
clear colors glme k normality whiskers dataset
%% Fig. S2b: Plot the mean-evoked traces from responders to each whisker in the D1 barrel at baseline
dataset = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 1,:);

Stm_C1 = []; Stm_D1 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(cell2mat(Stm_C1'), "omitnan"); std_C1 = std(cell2mat(Stm_C1'), "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(cell2mat(Stm_D1'), "omitnan"); std_D1 = std(cell2mat(Stm_D1'), "omitnan")/sqrt(width(Stm_D1));
clear Stm_C1 Stm_D1

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_C1', mean_D1']; stdErrors = [std_C1', std_D1'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-m','-b','-g','-y'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_C1,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of responsive cells in D1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means dataset greybar i style x
%% Fig. S2c: Plot the mean-evoked traces from non-responders to each whisker in the D1 barrel at baseline
dataset = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 0,:);

Stm_C1 = []; Stm_D1 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(cell2mat(Stm_C1'), "omitnan"); std_C1 = std(cell2mat(Stm_C1'), "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(cell2mat(Stm_D1'), "omitnan"); std_D1 = std(cell2mat(Stm_D1'), "omitnan")/sqrt(width(Stm_D1));
clear Stm_C1 Stm_D1

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_C1', mean_D1']; stdErrors = [std_C1', std_D1'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedNResp = figure;
axes_meanEvokedNResp = axes('Parent',figure_meanEvokedNResp);
hold(axes_meanEvokedNResp,'on');

style={'-m','-b','-g','-y'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_C1,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedNResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of non-responsive cells in D1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])


clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means  dataset greybar i style x
%% Fig. S2d: Plot the MI from Responders in the D1 barrel to all whiskers at baseline

% Choice of the dataset HERE
dataset = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 1,:);

Stm_C1 = []; Stm_D1 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(Stm_C1, "omitnan"); std_C1 = std(Stm_C1, "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(Stm_D1, "omitnan"); std_D1 = std(Stm_D1, "omitnan")/sqrt(width(Stm_D1));
clear Stm_C1 Stm_D1

xLabels = {'C1', 'D1'}; means = [mean_C1; mean_D1]; stdErrors = [std_C1; std_D1];
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];
% Plot the bar chart with error bars for the standard errors
figure_MIresp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MIresp.CData(k,:) = colors(k,:);
end
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:2
    scatter(k*(ones(length(dataset.MI(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.MI(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker stimulated'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of the responsive ROIs in D1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset.Whisker_stim = removecats(dataset.Whisker_stim, {'Stm_B2', 'Stm_E3'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means
clear adj_ci_cvrg adj_p h crit_p colors lme k normality whiskers dataset
%% Fig. S2e: Plot the MI from Non-Responders in the D1 barrel to all whiskers at baseline

% Choice of the dataset HERE
dataset = Data(Data.Barrel_imaged == "D1" & Data.Timepoint == "BL" & Data.ResponsivityStim == 0,:);

Stm_C1 = []; Stm_D1 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.Whisker_stim(row) == "Stm_C1" 
        Stm_C1 = [Stm_C1 dataset.MI(row)]; % change the variable
    elseif dataset.Whisker_stim(row) == "Stm_D1" 
        Stm_D1 = [Stm_D1 dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_C1 = mean(Stm_C1, "omitnan"); std_C1 = std(Stm_C1, "omitnan")/sqrt(width(Stm_C1)); % computation of the std error deviation
mean_D1 = mean(Stm_D1, "omitnan"); std_D1 = std(Stm_D1, "omitnan")/sqrt(width(Stm_D1));
clear Stm_C1 Stm_D1

xLabels = {'C1', 'D1'}; means = [mean_C1; mean_D1]; stdErrors = [std_C1; std_D1];
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];
% Plot the bar chart with error bars for the standard errors
figure_MIresp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MIresp.CData(k,:) = colors(k,:);
end
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'});
whiskers=unique(dataset.Whisker_stim);
for k = 1:2
    scatter(k*(ones(length(dataset.MI(dataset.Whisker_stim==whiskers(k))),1)),...
        dataset.MI(dataset.Whisker_stim==whiskers(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker stimulated'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of the non-responsive ROIs in D1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset.Whisker_stim = removecats(dataset.Whisker_stim, {'Stm_B2', 'Stm_E3'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means
clear adj_ci_cvrg adj_p h crit_p colors lme k normality whiskers dataset

%% Fig. S2f: Plot spatial distribution of responsive cells in D1 barrel at baseline

%Plot C1 Whisker
dataset = D1_BL(D1_BL.Whisker_stim == "Stm_C1",:); 

f_D1_stmC1_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color='m';
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('D1 Barrel - C1 Whisker Stimulation - Baseline');

%Plot D1 Whisker
dataset = D1_BL(D1_BL.Whisker_stim == "Stm_D1",:); 

f_D1_stmD1_BL_spatial=figure;
xlim([0 512])
ylim([0 512])
axis square
hold on

for i=1:size(dataset,1)
    center=fliplr(dataset.ROI_centers{i});
    radius=dataset.MI(i);
    if radius < 1
        radius = 1;
    else
    end
    
    if dataset.ResponsivityStim(i)==1
        color='b';
    else
        color='k';
    end
    drawcircle('Center',center,'Radius',radius,'Color',color,...
        'FaceAlpha',0.5,'StripeColor','none','HandleVisibility','off',...
        'InteractionsAllowed','none','LineWidth',0.001,'EdgeAlpha',0.5);
end
set(gca, 'YDir', 'reverse');
title('D1 Barrel - D1 Whisker Stimulation - Baseline');

clear center color dataset i radius
%% Fig. 2b: Plot selectivity index in C1 & D1 barrels at baseline

histogram(C1_BL_D46_all.SI_MI_BL, 20,"FaceColor", "m","FaceAlpha",0.5)
hold on
histogram(D1_BL_D46.SI_MI_BL, 20,"FaceColor", "b","FaceAlpha",0.5)
ylim([0 40])

C1_mean = mean(C1_BL_D46_all.SI_MI_BL); D1_mean = mean(D1_BL_D46.SI_MI_BL); yLimits = ylim;
C1_sem = std(C1_BL_D46_all.SI_MI_BL)/sqrt(length(C1_BL_D46_all.SI_MI_BL));
D1_sem = std(D1_BL_D46.SI_MI_BL)/sqrt(length(D1_BL_D46.SI_MI_BL));

line([C1_mean, C1_mean], yLimits, 'Color', 'm', 'LineStyle',':', 'LineWidth', 5);
line([D1_mean, D1_mean], yLimits, 'Color', 'b', 'LineStyle',':', 'LineWidth', 5); 
hold off

xlim([-1.1 1.1])

legend({'C1', 'D1', 'Mean C1', 'Mean D1'},'location','northwest')
title('SI in C1 & D1 barrel at baseline')

data_stats=[C1_BL_D46_all(:,{'MouseID','Barrel_imaged','Index','SI_MI_BL'});D1_BL_D46(:,{'MouseID','Barrel_imaged','Index','SI_MI_BL'})];
normality=lillietest(data_stats.SI_MI_BL) %1, data is not normal
lme = fitlme(data_stats,'SI_MI_BL ~ Barrel_imaged + (1| MouseID:Index)') %random intercept w/nested cells

clear lme lme2 normality data_stats
clear BL_mean D46 mean yLimits i row C1_mean C1_sem D1_mean D1_sem
%% Fig 2c: Plot changes in % Responsive cells to D1 whisker over time
%To generate figure S3b, remove the MI criteria in section 2, then run this
%section
mean_ByGroup=sortrows(mean_ByGroup,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_PrctRespChronic = figure;
axes_PrctRespChronic = axes('Parent',figure_PrctRespChronic);
hold(axes_PrctRespChronic,'on');

errorbar(Timepoints,...
    (mean_ByGroup.PercentResponsive(mean_ByGroup.Barrel_imaged=='D1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    (mean_ByGroup.PercentResponsive_SEM(mean_ByGroup.Barrel_imaged=='D1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    'DisplayName','D2 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (mean_ByGroup.PercentResponsive(mean_ByGroup.Barrel_imaged=='C1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    (mean_ByGroup.PercentResponsive_SEM(mean_ByGroup.Barrel_imaged=='C1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    'DisplayName','C1 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_PrctRespChronic,'FontSize',12,'YTick',[0 20 40 60 80 100],'ylim',[0 100],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_PrctRespChronic, '% Responsive Cells')
title(axes_PrctRespChronic,'% D1 Responsive cells in C1 and D1 barrels')

% Add x-axis labels and title
xticks(Timepoints);
xLabels = {'0', '3', '6', '12', '18', '46'};
xlabel('Time (days)');
xticklabels(xLabels);

%Get statistics
dataset = mean_ByMouse(mean_ByMouse.Whisker_stim=='Stm_D1' & mean_ByMouse.Barrel_imaged=='D1',:);
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Timepoint + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))
aov = anova(glme)
clear adj_p adj_ci_cvrg crit_p normality aov glme h

dataset = mean_ByMouse(mean_ByMouse.Whisker_stim=='Stm_D1' & mean_ByMouse.Barrel_imaged=='C1',:);
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Timepoint + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))
aov = anova(glme)

clear adj_p adj_ci_cvrg crit_p normality aov glme h
clear  xLabels dataset Timepoints

%% Fig 2d: Plot MI for D1-Responsive Cells in C1 and D1 barrels over time
% To generate figure S3c, remove the MI criteria in Section 2, then run
% this section
mean_Resp=sortrows(mean_Resp,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_MI_r = figure;
axes_MI_r = axes('Parent',figure_MI_r);
hold(axes_MI_r,'on');

errorbar(Timepoints,...
    (mean_Resp.MI(mean_Resp.Barrel_imaged=='D1' & mean_Resp.Whisker_stim=='Stm_D1')),...
    (mean_Resp.MI_SEM(mean_Resp.Barrel_imaged=='D1' & mean_Resp.Whisker_stim=='Stm_D1')),...
    'DisplayName','Sham D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (mean_Resp.MI(mean_Resp.Barrel_imaged=='C1' & mean_Resp.Whisker_stim=='Stm_D1')),...
    (mean_Resp.MI_SEM(mean_Resp.Barrel_imaged=='C1' & mean_Resp.Whisker_stim=='Stm_D1')),...
    'DisplayName','Stroke D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_MI_r,'FontSize',12,'YTick',[0 2 4 6 8 10],'ylim',[0 10],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_MI_r, 'MI (Z-score)')
title(axes_MI_r,'MI - D1-Responsive PV cells')

%Get statistics for D1 barrel
dataset = Data_Resp(Data_Resp.Whisker_stim=='Stm_D1' & Data_Resp.Barrel_imaged=='D1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal
% lme = fitlme(dataset,'MI ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))
clear h adj_p adj_ci_cvrg crit_p normality

%Get statistics for C1 barrel
dataset = Data_Resp(Data_Resp.Whisker_stim=='Stm_D1' & Data_Resp.Barrel_imaged=='C1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal
% lme = fitlme(dataset,'MI ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))
clear h adj_p adj_ci_cvrg crit_p normality

clear  xLabels aov dataset lme2 Timepoints

%% Fig. 2e: Plot change in responsivity of cells in D1 barrel to D1 stim - BL vs D12
% Carlos Borau (2024). Sankey flow chart (https://www.mathworks.com/matlabcentral/fileexchange/101516-sankey-flow-chart), MATLAB Central File Exchange. Retrieved January 3, 2024.
% Colormap: can be the name of matlab colormaps or a matrix of (N x 3).
%   Important: N must be the max number of categories in a layer 
%   multiplied by the number of layers. 
%   In the example of Option 1, N should be 4 * 4 = 16;
%   In the example of Option 2, N should be 5 * 3 = 15;
Sandata=D1_Stm_D1(:,[12,15]);         
Sandata.Responsiveness_BL_D1=categorical(Sandata.Responsiveness_BL_D1);
Sandata.Responsiveness_D12_D1=categorical(Sandata.Responsiveness_D12_D1);
Sandata.Responsiveness_BL_D1=renamecats(Sandata.Responsiveness_BL_D1,{'Non-responsive','Responsive'});
Sandata.Responsiveness_D12_D1=renamecats(Sandata.Responsiveness_D12_D1,{'Non-responsive','Responsive'});
Sandata=sortrows(Sandata);
options.color_map = [0.3 0.3 0.3; 1 0 0; 0.3 0.3 0.3; 1 0 0];    
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = false;          % show percentage over the blocks
options.text_color = [1 1 1];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;
plotSankeyFlowChart(Sandata,options);
title('Change in Responsivity of Cells in D1 barrel to D1 - BL vs D12');
clear Sandata options
%% Fig. 2f: Plot change in responsivity of cells in C1 barrel to D1 stim - BL vs D12
% Carlos Borau (2024). Sankey flow chart (https://www.mathworks.com/matlabcentral/fileexchange/101516-sankey-flow-chart), MATLAB Central File Exchange. Retrieved January 3, 2024.
% Colormap: can be the name of matlab colormaps or a matrix of (N x 3).
%   Important: N must be the max number of categories in a layer 
%   multiplied by the number of layers. 
%   In the example of Option 1, N should be 4 * 4 = 16;
%   In the example of Option 2, N should be 5 * 3 = 15;

Sandata=C1_Stm_D1(:,[12,15]);         
Sandata.Responsiveness_BL_D1=categorical(Sandata.Responsiveness_BL_D1);
Sandata.Responsiveness_D12_D1=categorical(Sandata.Responsiveness_D12_D1);
Sandata.Responsiveness_BL_D1=renamecats(Sandata.Responsiveness_BL_D1,{'Non-responsive','Responsive'});
Sandata.Responsiveness_D12_D1=renamecats(Sandata.Responsiveness_D12_D1,{'Non-responsive','Responsive'});
Sandata=sortrows(Sandata);
options.color_map = [0.3 0.3 0.3; 1 0 0; 0.3 0.3 0.3; 1 0 0];      
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = false;          % show percentage over the blocks
options.text_color = [1 1 1];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;
plotSankeyFlowChart(Sandata,options);
title('Change in Responsivity of Cells in C1 barrel to D1 - BL vs D12');
clear Sandata options

%% Fig. 2g: Plot the change in MI in D1 barrel to D1 stim from baseline to d12 by switch category
dataset= D1_Stm_D1;
dataset.dMI_BLtoD12=dataset.MI_D12_D1-dataset.MI_BL_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12)
NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.dMI_BLtoD12(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));
figure
xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.dMI_BLtoD12(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.dMI_BLtoD12(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Change in MI of Cells in D1 barrel to D1 stim BL to D12'); % Choice of the name of the plot HERE
ylim([-30 50]) % add a limit here if needed to scale different plots together reliability 

normality = lillietest(dataset.dMI_BLtoD12) %1, not normal
[p,tbl,stats] = kruskalwallis(dataset.dMI_BLtoD12,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    k colors means p stats tbl normality c

%% Fig. 2h: Plot the change in MI in C1 barrel to D1 stim from baseline to d12 by switch category
dataset= C1_Stm_D1;
dataset.dMI_BLtoD12=dataset.MI_D12_D1-dataset.MI_BL_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12)
NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.dMI_BLtoD12(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.dMI_BLtoD12(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));
figure
xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.dMI_BLtoD12(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.dMI_BLtoD12(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Change in MI of Cells in C1 barrel to D1 stim BL to D12'); % Choice of the name of the plot HERE
ylim([-30 50]) % add a limit here if needed to scale different plots together reliability 

normality = lillietest(dataset.dMI_BLtoD12) %1, not normal
[p,tbl,stats] = kruskalwallis(dataset.dMI_BLtoD12,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    k colors means p stats tbl normality c

%% Fig. 2i: Plot MI of cells in D1 barrel to D1 stim @ D12 by switch category
dataset = D1_Stm_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12)
NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI_D12_D1(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));
figure
xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_D12_byswitch=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_D12_byswitch.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI_D12_D1(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI_D12_D1(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in D1 barrel to D1 stim at D12'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability 

normality = lillietest(dataset.MI_D12_D1) %1, not normal
[p,tbl,stats] = kruskalwallis(dataset.MI_D12_D1,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means p stats tbl normality c
%% Fig. 2j: Plot MI of cells in C1 barrel to D1 stim @ D12 by switch category
dataset = C1_Stm_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12)
NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI_D12_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI_D12_D1(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));
figure
xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_D12_byswitch_C1=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_D12_byswitch_C1.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI_D12_D1(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI_D12_D1(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in C1 barrel to D1 stim at D12'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability 

normality = lillietest(dataset.MI_D12_D1) %1, not normal
[p,tbl,stats] = kruskalwallis(dataset.MI_D12_D1,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means p stats tbl normality c
%% Fig. S3b: Plot %D1 responsive cells in C1 and D1 barrels over time without using MI threshold
%First repeat sections 1 & 2, removing the MI theshold criteria.
%Then repeat section 3.
%Then run the section to plot Figure 2c

%% Fig. S3c: Plot Changes in MI to D1 whisker in C1 and D1 barrels over time without using MI threshold
%First repeat sections 1 & 2, removing the MI theshold criteria.
%Then repeat section 3.
%Then run the section to plot Figure 2d
%% Fig. S3d: Plot MI of cells in D1 barrel to D1 stim over time
dataset = D1_Stm_D1;

mean_BL = mean(D1_Stm_D1.MI_BL_D1, "omitnan"); std_BL = std(D1_Stm_D1.MI_BL_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_BL_D1)); % computation of the std error deviation
mean_d3 = mean(D1_Stm_D1.MI_D3_D1, "omitnan"); std_d3 = std(D1_Stm_D1.MI_D3_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_D3_D1));
mean_d6 = mean(D1_Stm_D1.MI_D6_D1, "omitnan"); std_d6 = std(D1_Stm_D1.MI_D6_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_D6_D1));
mean_d12 = mean(D1_Stm_D1.MI_D12_D1, "omitnan"); std_d12 = std(D1_Stm_D1.MI_D12_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_D12_D1));
mean_d18 = mean(D1_Stm_D1.MI_D18_D1, "omitnan"); std_d18 = std(D1_Stm_D1.MI_D18_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_D18_D1));
mean_d46 = mean(D1_Stm_D1.MI_D46_D1, "omitnan"); std_d46 = std(D1_Stm_D1.MI_D46_D1, "omitnan")/sqrt(length(D1_Stm_D1.MI_D46_D1));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];

% Plot the bar chart with error bars for the standard errors
figure_MI_D1_allcells=bar(means,'FaceColor',[0 0 1],'FaceAlpha',0.5);
hold on

for k = 1:size(means,1)
    scatter(k*(ones(height(dataset),1)),dataset{:,k+5},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','b','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in D1 barrel to D1 stim over time'); % Choice of the name of the plot HERE
ylim([-5 80]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for D1 barrel
dataset = Data(Data.Whisker_stim=='Stm_D1' & Data.Barrel_imaged=='D1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal
% lme = fitlme(dataset,'MI ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    k means aov lme2 stdErrors

%% Fig. S3e: Plot MI of cells in C1 barrel to D1 stim over time
dataset = C1_Stm_D1;

mean_BL = mean(C1_Stm_D1.MI_BL_D1, "omitnan"); std_BL = std(C1_Stm_D1.MI_BL_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_BL_D1)); % computation of the std error deviation
mean_d3 = mean(C1_Stm_D1.MI_D3_D1, "omitnan"); std_d3 = std(C1_Stm_D1.MI_D3_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_D3_D1));
mean_d6 = mean(C1_Stm_D1.MI_D6_D1, "omitnan"); std_d6 = std(C1_Stm_D1.MI_D6_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_D6_D1));
mean_d12 = mean(C1_Stm_D1.MI_D12_D1, "omitnan"); std_d12 = std(C1_Stm_D1.MI_D12_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_D12_D1));
mean_d18 = mean(C1_Stm_D1.MI_D18_D1, "omitnan"); std_d18 = std(C1_Stm_D1.MI_D18_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_D18_D1));
mean_d46 = mean(C1_Stm_D1.MI_D46_D1, "omitnan"); std_d46 = std(C1_Stm_D1.MI_D46_D1, "omitnan")/sqrt(length(C1_Stm_D1.MI_D46_D1));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];

% Plot the bar chart with error bars for the standard errors
figure_MI_C1_allcells=bar(means,'FaceColor',[1 0 1],'FaceAlpha',0.5);
hold on

for k = 1:size(means,1)
    scatter(k*(ones(height(dataset),1)),dataset{:,k+5},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','m','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in C1 barrel to D1 stim over time'); % Choice of the name of the plot HERE
ylim([-5 80]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for C1 barrel
dataset = Data(Data.Whisker_stim=='Stm_D1' & Data.Barrel_imaged=='C1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal
% lme = fitlme(dataset,'MI ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    k means aov lme2 stdErrors
%% Fig. S3f: Calculate the mean %PV cells co-modulated (MI>1) during whisker stimulation
%First calculate the percentage of whisker responsive cells within a mouse
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMovie(:,{'MouseID','Timepoint','Barrel_imaged','Whisker_stim'});
[G]=findgroups(T1);

mean_ByMouse.stim_PrctCellMod=splitapply(omean,mean_ByMovie.stim_PrctCellMod,G);

%Then summarize by group to get mean with n=# of mice
T1=mean_ByMouse(:,{'Timepoint','Barrel_imaged','Whisker_stim'});
[G]=findgroups(T1);

mean_ByGroup.stim_PrctCellMod=splitapply(omean,mean_ByMouse.stim_PrctCellMod,G);
mean_ByGroup.stim_PrctCellMod_SEM=splitapply(sem,mean_ByMouse.stim_PrctCellMod,G);

clear PrctCellMod_ByMouse PrctCellMod_ByGroup

%Plot the graph
mean_ByGroup=sortrows(mean_ByGroup,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_PrctModChronic = figure;
axes_PrctModChronic = axes('Parent',figure_PrctModChronic);
hold(axes_PrctModChronic,'on');

errorbar(Timepoints,...
    (mean_ByGroup.stim_PrctCellMod(mean_ByGroup.Barrel_imaged=='D1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    (mean_ByGroup.stim_PrctCellMod_SEM(mean_ByGroup.Barrel_imaged=='D1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    'DisplayName','D2 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (mean_ByGroup.stim_PrctCellMod(mean_ByGroup.Barrel_imaged=='C1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    (mean_ByGroup.stim_PrctCellMod_SEM(mean_ByGroup.Barrel_imaged=='C1' & mean_ByGroup.Whisker_stim=='Stm_D1')),...
    'DisplayName','C1 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_PrctModChronic,'FontSize',12,'YTick',[0 .2 .4 .6 .8 1],'ylim',[0 1],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_PrctModChronic, '% Modulated Cells/Stimulus')
title(axes_PrctModChronic,'% Modulated PV cells in C1 and D1 barrels per D1 whisker stimulus')

% Add x-axis labels and title
xticks(Timepoints);
xLabels = {'0', '3', '6', '12', '18', '46'};
xlabel('Time (days)');
xticklabels(xLabels);

%Get statistics
dataset = mean_ByMouse(mean_ByMouse.Whisker_stim=='Stm_D1' & mean_ByMouse.Barrel_imaged=='D1',:);
normality = lillietest(dataset.stim_PrctCellMod) %1, not normal
lme = fitlme(dataset,'stim_PrctCellMod ~ Timepoint + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov = anova(lme)
clear adj_p adj_ci_cvrg crit_p normality aov lme h

dataset = mean_ByMouse(mean_ByMouse.Whisker_stim=='Stm_D1' & mean_ByMouse.Barrel_imaged=='C1',:);
normality = lillietest(dataset.stim_PrctCellMod) 
lme = fitlme(dataset,'stim_PrctCellMod ~ Timepoint + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov = anova(lme)

clear adj_p adj_ci_cvrg crit_p normality aov lme h
clear  xLabels dataset Timepoints omean sem G T1
%% Fig. S3g: Plot Xcorr of cells in D1 barrel to D1 stim over time
dataset = Data_ByCell_wCorr(Data_ByCell_wCorr.Barrel_imaged=='D1' &...
    Data_ByCell_wCorr.Whisker_stim=='Stm_D1',:);

dataset=dataset(~isnan(dataset.Rreal_avg),:); %remove cells where xcorr could not be correlated because there was only 1 cell

%Get means for true cross correlation
BL = table2array(dataset(dataset.Timepoint=='BL','Rreal_avg')); 
d3 = table2array(dataset(dataset.Timepoint=='Day3','Rreal_avg')); 
d6 = table2array(dataset(dataset.Timepoint=='Day6','Rreal_avg'));  
d12 = table2array(dataset(dataset.Timepoint=='Day12','Rreal_avg')); 
d18 = table2array(dataset(dataset.Timepoint=='Day18','Rreal_avg')); 
d46 = table2array(dataset(dataset.Timepoint=='Day46','Rreal_avg'));
alltimes=[BL,d3,d6,d12,d18,d46];

mean_BL = mean(BL, "omitnan"); std_BL = std(BL, "omitnan")/sqrt(length(BL)); % computation of the std error deviation
mean_d3 = mean(d3, "omitnan"); std_d3 = std(d3, "omitnan")/sqrt(length(d3));
mean_d6 = mean(d6, "omitnan"); std_d6 = std(d6, "omitnan")/sqrt(length(d6));
mean_d12 = mean(d12, "omitnan"); std_d12 = std(d12, "omitnan")/sqrt(length(d12));
mean_d18 = mean(d18, "omitnan"); std_d18 = std(d18, "omitnan")/sqrt(length(d18));
mean_d46 = mean(d46, "omitnan"); std_d46 = std(d46, "omitnan")/sqrt(length(d46));

%Get the means for shuffled cross correlation
BLshuff = table2array(dataset(dataset.Timepoint=='BL','Rshuff_avg')); 
d3shuff = table2array(dataset(dataset.Timepoint=='Day3','Rshuff_avg')); 
d6shuff = table2array(dataset(dataset.Timepoint=='Day6','Rshuff_avg'));  
d12shuff = table2array(dataset(dataset.Timepoint=='Day12','Rshuff_avg')); 
d18shuff = table2array(dataset(dataset.Timepoint=='Day18','Rshuff_avg')); 
d46shuff = table2array(dataset(dataset.Timepoint=='Day46','Rshuff_avg'));
alltimesshuff=[BL,d3,d6,d12,d18,d46];

mean_BLshuff = mean(BLshuff, "omitnan"); std_BLshuff = std(BLshuff, "omitnan")/sqrt(length(BLshuff)); % computation of the std error deviation
mean_d3shuff = mean(d3shuff, "omitnan"); std_d3shuff = std(d3shuff, "omitnan")/sqrt(length(d3shuff));
mean_d6shuff = mean(d6shuff, "omitnan"); std_d6shuff = std(d6shuff, "omitnan")/sqrt(length(d6shuff));
mean_d12shuff = mean(d12shuff, "omitnan"); std_d12shuff = std(d12shuff, "omitnan")/sqrt(length(d12shuff));
mean_d18shuff = mean(d18shuff, "omitnan"); std_d18shuff = std(d18shuff, "omitnan")/sqrt(length(d18shuff));
mean_d46shuff = mean(d46shuff, "omitnan"); std_d46shuff = std(d46shuff, "omitnan")/sqrt(length(d46shuff));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];
means_shuff = [mean_BLshuff; mean_d3shuff; mean_d6shuff; mean_d12shuff; mean_d18shuff; mean_d46shuff]; 
stdErrors_shuff = [std_BLshuff; std_d3shuff; std_d6shuff; std_d12shuff; std_d18shuff; std_d46shuff];

% Plot the bar chart with error bars for the standard errors
figure_xorr_D1_allcells=bar(means,'FaceColor',[0 0 1],'FaceAlpha',0.5);
hold on
bar(means_shuff,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.5);

for k = 1:size(means,1)
    scatter(k*(ones(height(alltimes),1)),alltimes(:,k),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','b','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
errorbar(means_shuff, stdErrors_shuff, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('Xcorr - (R-value)'); % Choice of the name of the variable plot HERE
title('Cross Correlation of Cells in D1 barrel to D1 stim over time'); % Choice of the name of the plot HERE
ylim([-0.4 1]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for D1 barrel
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.Rreal_avg) %1=not normal
% lme = fitlme(dataset,'Rreal_avg ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'Rreal_avg ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality aov lme2
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    BL d3 d6 d12 d18 d46 BLshuff d3shuff d6shuff d12shuff d18shuff d46shuff...
    k means mean_BLshuff mean_d3shuff mean_d6shuff mean_d12shuff mean_d18shuff mean_d46shuff...
    std_BLshuff std_d3shuff std_d6shuff std_d12shuff std_d18shuff std_d46shuff...
    means_shuff stdErrors_shuff alltimes alltimesshuff stdErrors

%% Fig. S3h: Plot Xcorr of cells in C1 barrel to D1 stim over time
dataset = Data_ByCell_wCorr(Data_ByCell_wCorr.Barrel_imaged=='C1' &...
    Data_ByCell_wCorr.Whisker_stim=='Stm_D1',:);

dataset=dataset(~isnan(dataset.Rreal_avg),:); %remove cells where xcorr could not be correlated because there was only 1 cell

%Get means for true cross correlation
BL = table2array(dataset(dataset.Timepoint=='BL','Rreal_avg')); 
d3 = table2array(dataset(dataset.Timepoint=='Day3','Rreal_avg')); 
d6 = table2array(dataset(dataset.Timepoint=='Day6','Rreal_avg'));  
d12 = table2array(dataset(dataset.Timepoint=='Day12','Rreal_avg')); 
d18 = table2array(dataset(dataset.Timepoint=='Day18','Rreal_avg')); 
d46 = table2array(dataset(dataset.Timepoint=='Day46','Rreal_avg'));
alltimes=[BL,d3,d6,d12,d18,d46];

mean_BL = mean(BL, "omitnan"); std_BL = std(BL, "omitnan")/sqrt(length(BL)); % computation of the std error deviation
mean_d3 = mean(d3, "omitnan"); std_d3 = std(d3, "omitnan")/sqrt(length(d3));
mean_d6 = mean(d6, "omitnan"); std_d6 = std(d6, "omitnan")/sqrt(length(d6));
mean_d12 = mean(d12, "omitnan"); std_d12 = std(d12, "omitnan")/sqrt(length(d12));
mean_d18 = mean(d18, "omitnan"); std_d18 = std(d18, "omitnan")/sqrt(length(d18));
mean_d46 = mean(d46, "omitnan"); std_d46 = std(d46, "omitnan")/sqrt(length(d46));

%Get the means for shuffled cross correlation
BLshuff = table2array(dataset(dataset.Timepoint=='BL','Rshuff_avg')); 
d3shuff = table2array(dataset(dataset.Timepoint=='Day3','Rshuff_avg')); 
d6shuff = table2array(dataset(dataset.Timepoint=='Day6','Rshuff_avg'));  
d12shuff = table2array(dataset(dataset.Timepoint=='Day12','Rshuff_avg')); 
d18shuff = table2array(dataset(dataset.Timepoint=='Day18','Rshuff_avg')); 
d46shuff = table2array(dataset(dataset.Timepoint=='Day46','Rshuff_avg'));
alltimesshuff=[BL,d3,d6,d12,d18,d46];

mean_BLshuff = mean(BLshuff, "omitnan"); std_BLshuff = std(BLshuff, "omitnan")/sqrt(length(BLshuff)); % computation of the std error deviation
mean_d3shuff = mean(d3shuff, "omitnan"); std_d3shuff = std(d3shuff, "omitnan")/sqrt(length(d3shuff));
mean_d6shuff = mean(d6shuff, "omitnan"); std_d6shuff = std(d6shuff, "omitnan")/sqrt(length(d6shuff));
mean_d12shuff = mean(d12shuff, "omitnan"); std_d12shuff = std(d12shuff, "omitnan")/sqrt(length(d12shuff));
mean_d18shuff = mean(d18shuff, "omitnan"); std_d18shuff = std(d18shuff, "omitnan")/sqrt(length(d18shuff));
mean_d46shuff = mean(d46shuff, "omitnan"); std_d46shuff = std(d46shuff, "omitnan")/sqrt(length(d46shuff));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];
means_shuff = [mean_BLshuff; mean_d3shuff; mean_d6shuff; mean_d12shuff; mean_d18shuff; mean_d46shuff]; 
stdErrors_shuff = [std_BLshuff; std_d3shuff; std_d6shuff; std_d12shuff; std_d18shuff; std_d46shuff];

% Plot the bar chart with error bars for the standard errors
figure_xcorr_C1_allcells=bar(means,'FaceColor',[1 0 1],'FaceAlpha',0.5);
hold on
bar(means_shuff,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.5);

for k = 1:size(means,1)
    scatter(k*(ones(height(alltimes),1)),alltimes(:,k),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','m','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
errorbar(means_shuff, stdErrors_shuff, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('Xcorr - (R-value)'); % Choice of the name of the variable plot HERE
title('Cross Correlation of Cells in C1 barrel to D1 stim over time'); % Choice of the name of the plot HERE
ylim([-0.4 1]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for C1 barrel
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.Rreal_avg) %1=not normal
% lme = fitlme(dataset,'Rreal_avg ~ Timepoint + (1 | MouseID:Index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'Rreal_avg ~ Timepoint + (Timepoint | MouseID:Index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality aov lme2
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    BL d3 d6 d12 d18 d46 BLshuff d3shuff d6shuff d12shuff d18shuff d46shuff...
    k means mean_BLshuff mean_d3shuff mean_d6shuff mean_d12shuff mean_d18shuff mean_d46shuff...
    std_BLshuff std_d3shuff std_d6shuff std_d12shuff std_d18shuff std_d46shuff...
    means_shuff stdErrors_shuff alltimes alltimesshuff stdErrors
%% Fig. S4a: No Trim Control - Plot changes in % Responsive cells to D1 whisker over time
mean_ByGroup_NoTrim=sortrows(mean_ByGroup_NoTrim,{'Timepoint'});
Timepoints=[0,18,46];
figure_PrctRespChronic = figure;
axes_PrctRespChronic = axes('Parent',figure_PrctRespChronic);
hold(axes_PrctRespChronic,'on');

errorbar(Timepoints,...
    (mean_ByGroup_NoTrim.PercentResponsive(mean_ByGroup_NoTrim.Barrel_imaged=='D1' & mean_ByGroup_NoTrim.Whisker_stim=='D1')),...
    (mean_ByGroup_NoTrim.PercentResponsive_SEM(mean_ByGroup_NoTrim.Barrel_imaged=='D1' & mean_ByGroup_NoTrim.Whisker_stim=='D1')),...
    'DisplayName','D1 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (mean_ByGroup_NoTrim.PercentResponsive(mean_ByGroup_NoTrim.Barrel_imaged=='C1' & mean_ByGroup_NoTrim.Whisker_stim=='D1')),...
    (mean_ByGroup_NoTrim.PercentResponsive_SEM(mean_ByGroup_NoTrim.Barrel_imaged=='C1' & mean_ByGroup_NoTrim.Whisker_stim=='D1')),...
    'DisplayName','C1 Barrel','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_PrctRespChronic,'FontSize',12,'YTick',[0 20 40 60 80 100],'ylim',[0 100],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_PrctRespChronic, '% Responsive Cells')
title(axes_PrctRespChronic,'% D1 Responsive cells in C1 and D1 barrels')

% Add x-axis labels and title
xticks(Timepoints);
xLabels = {'0', '18', '46'};
xlabel('Time (days)');
xticklabels(xLabels);

%Get statistics
dataset = mean_ByMouse_NoTrim(mean_ByMouse_NoTrim.Whisker_stim=='D1' & mean_ByMouse_NoTrim.Barrel_imaged=='D1',:);
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Timepoint + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))
aov = anova(glme)
clear adj_p adj_ci_cvrg crit_p normality

dataset = mean_ByMouse_NoTrim(mean_ByMouse_NoTrim.Whisker_stim=='D1' & mean_ByMouse_NoTrim.Barrel_imaged=='C1',:);
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Timepoint + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))
aov = anova(glme)
clear adj_p adj_ci_cvrg crit_p normality

clear  xLabels aov glme h Timepoints

%% Fig. S4b: No Trim Control - Plot MI for D1-Responsive Cells in C1 and D1 barrels over time
mean_Resp_NoTrim=sortrows(mean_Resp_NoTrim,{'Timepoint'});
Timepoints=[0,18,46];
figure_MI_r = figure;
axes_MI_r = axes('Parent',figure_MI_r);
hold(axes_MI_r,'on');

errorbar(Timepoints,...
    (mean_Resp_NoTrim.MI(mean_Resp_NoTrim.Barrel_imaged=='D1' & mean_Resp_NoTrim.Whisker_stim=='D1')),...
    (mean_Resp_NoTrim.MI_SEM(mean_Resp_NoTrim.Barrel_imaged=='D1' & mean_Resp_NoTrim.Whisker_stim=='D1')),...
    'DisplayName','Sham D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (mean_Resp_NoTrim.MI(mean_Resp_NoTrim.Barrel_imaged=='C1' & mean_Resp_NoTrim.Whisker_stim=='D1')),...
    (mean_Resp_NoTrim.MI_SEM(mean_Resp_NoTrim.Barrel_imaged=='C1' & mean_Resp_NoTrim.Whisker_stim=='D1')),...
    'DisplayName','Stroke D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_MI_r,'FontSize',12,'YTick',[0 4 8 12],'ylim',[0 12],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_MI_r, 'MI (Z-score)')
title(axes_MI_r,'MI - D1-Responsive PV cells')
% Add x-axis labels and title
xticks(Timepoints);
xLabels = {'0', '18', '46'};
xlabel('Time (days)');
xticklabels(xLabels);

%Get statistics for D1 barrel
dataset = Data_Resp_NoTrim(Data_Resp_NoTrim.Whisker_stim=='D1' & Data_Resp_NoTrim.Barrel_imaged=='D1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal

lme = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:index)') %random slope & intercept
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality

%Get statistics for C1 barrel
dataset = Data_Resp_NoTrim(Data_Resp_NoTrim.Whisker_stim=='D1' & Data_Resp_NoTrim.Barrel_imaged=='C1',:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim); 
normality = lillietest(dataset.MI) %1=not normal

lme = fitlme(dataset,'MI ~ Timepoint + (Timepoint | MouseID:index)') %random slope & intercept
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality

clear  xLabels aov lme Timepoints
%% Fig. S5a: Plot AUC of spontaneous activity for All Cells in C1 and D1 barrels over time
spont_mean_allcells=sortrows(spont_mean_allcells,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_AUC = figure;
axes_AUC = axes('Parent',figure_AUC);
hold(axes_AUC,'on');

errorbar(Timepoints,...
    (spont_mean_allcells.AUC(spont_mean_allcells.Barrel_imaged=='D1')),...
    (spont_mean_allcells.AUC_SEM(spont_mean_allcells.Barrel_imaged=='D1')),...
    'DisplayName','Sham D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (spont_mean_allcells.AUC(spont_mean_allcells.Barrel_imaged=='C1')),...
    (spont_mean_allcells.AUC_SEM(spont_mean_allcells.Barrel_imaged=='C1')),...
    'DisplayName','Stroke D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_AUC,'FontSize',12,'YTick',[0 100 200 300 400 500],'ylim',[0 500],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_AUC, 'AUC (Z-score * sec)')
title(axes_AUC,'AUC PV cell Spontaneous Activity')

%Get statistics for D1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='D1',:);
normality = lillietest(dataset.AUC) %1=not normal

lme = fitlme(dataset,'AUC ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality

%Get statistics for C1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='C1',:);
normality = lillietest(dataset.AUC) %1=not normal

lme = fitlme(dataset,'AUC ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality aov row Timepoints lme
%% Fig. S5b: Plot Peak Frequency of spontaneous activity for All Cells in C1 and D1 barrels over time
spont_mean_allcells=sortrows(spont_mean_allcells,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_peakfreq = figure;
axes_peakfreq = axes('Parent',figure_peakfreq);
hold(axes_peakfreq,'on');

errorbar(Timepoints,...
    (spont_mean_allcells.PeakFreq(spont_mean_allcells.Barrel_imaged=='D1')),...
    (spont_mean_allcells.PeakFreq_SEM(spont_mean_allcells.Barrel_imaged=='D1')),...
    'DisplayName','Sham D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (spont_mean_allcells.PeakFreq(spont_mean_allcells.Barrel_imaged=='C1')),...
    (spont_mean_allcells.PeakFreq_SEM(spont_mean_allcells.Barrel_imaged=='C1')),...
    'DisplayName','Stroke D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_peakfreq,'FontSize',12,'YTick',[0 0.1 0.2 0.3 0.4],'ylim',[0 0.4],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_peakfreq, 'Calcium Transient Frequency')
title(axes_peakfreq,'Mean Frequency of PV cell Spontaneous Activity')

%Get statistics for D1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='D1',:);
normality = lillietest(dataset.PeakFreq) %1=not normal

lme = fitlme(dataset,'PeakFreq ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality

%Get statistics for C1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='C1',:);
normality = lillietest(dataset.PeakFreq) %1=not normal

lme = fitlme(dataset,'PeakFreq ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality aov lme Timepoints

%% Fig. S5c: Plot Mean Peak Amplitude of spontaneous activity for All Cells in C1 and D1 barrels over time
spont_mean_allcells=sortrows(spont_mean_allcells,{'Timepoint'});
Timepoints=[0,3,6,12,18,46];
figure_peakamp = figure;
axes_peakamp = axes('Parent',figure_peakamp);
hold(axes_peakamp,'on');

errorbar(Timepoints,...
    (spont_mean_allcells.PeaksAmp_mean(spont_mean_allcells.Barrel_imaged=='D1')),...
    (spont_mean_allcells.PeaksAmp_mean_SEM(spont_mean_allcells.Barrel_imaged=='D1')),...
    'DisplayName','Sham D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','b',...
'Color','b')

errorbar(Timepoints,...
    (spont_mean_allcells.PeaksAmp_mean(spont_mean_allcells.Barrel_imaged=='C1')),...
    (spont_mean_allcells.PeaksAmp_mean_SEM(spont_mean_allcells.Barrel_imaged=='C1')),...
    'DisplayName','Stroke D2 - Stim D2','MarkerSize',10,...
'MarkerEdgeColor','none','Marker','o','LineStyle','-','LineWidth',3,'MarkerFaceColor','m',...
'Color','m')

set(axes_peakamp,'FontSize',12,'YTick',[0 2 4 6 8 10 12],'ylim',[0 12],'XTick',Timepoints,'xlim',[-5 51]);
% legend('Location','northeastoutside')
ylabel(axes_peakamp, 'Calcium Transient Amplitude (Z-score)')
title(axes_peakamp,'Mean Amplitude of PV cell Spontaneous Activity')

%Get statistics for D1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='D1',:);
normality = lillietest(dataset.PeaksAmp_mean) %1=not normal

lme = fitlme(dataset,'PeaksAmp_mean ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality

%Get statistics for C1 barrel
dataset = Spont_Data_ByCell(Spont_Data_ByCell.Barrel_imaged=='C1',:);
normality = lillietest(dataset.PeaksAmp_mean) %1=not normal

lme = fitlme(dataset,'PeaksAmp_mean ~ Timepoint + (index | MouseID)')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))
aov=anova(lme)
clear h adj_p adj_ci_cvrg crit_p normality aov lme Timepoints
%% Fig. S5d: Plot Xcorr of spontaneous activity of cells in D1 barrel over time
dataset = Spont_Data_ByCell_wCorr(Spont_Data_ByCell_wCorr.Barrel_imaged=='D1',:);

dataset=dataset(~isnan(dataset.Rreal_avg),:); %remove cells where xcorr could not be correlated because there was only 1 cell

%Get means for true cross correlation
BL = table2array(dataset(dataset.Timepoint=='0','Rreal_avg')); 
d3 = table2array(dataset(dataset.Timepoint=='3','Rreal_avg')); 
d6 = table2array(dataset(dataset.Timepoint=='6','Rreal_avg'));  
d12 = table2array(dataset(dataset.Timepoint=='12','Rreal_avg')); 
d18 = table2array(dataset(dataset.Timepoint=='18','Rreal_avg')); 
d46 = table2array(dataset(dataset.Timepoint=='46','Rreal_avg'));
alltimes=[BL,d3,d6,d12,d18,d46];

mean_BL = mean(BL, "omitnan"); std_BL = std(BL, "omitnan")/sqrt(length(BL)); % computation of the std error deviation
mean_d3 = mean(d3, "omitnan"); std_d3 = std(d3, "omitnan")/sqrt(length(d3));
mean_d6 = mean(d6, "omitnan"); std_d6 = std(d6, "omitnan")/sqrt(length(d6));
mean_d12 = mean(d12, "omitnan"); std_d12 = std(d12, "omitnan")/sqrt(length(d12));
mean_d18 = mean(d18, "omitnan"); std_d18 = std(d18, "omitnan")/sqrt(length(d18));
mean_d46 = mean(d46, "omitnan"); std_d46 = std(d46, "omitnan")/sqrt(length(d46));

%Get the means for shuffled cross correlation
BLshuff = table2array(dataset(dataset.Timepoint=='0','Rshuff_avg')); 
d3shuff = table2array(dataset(dataset.Timepoint=='3','Rshuff_avg')); 
d6shuff = table2array(dataset(dataset.Timepoint=='6','Rshuff_avg'));  
d12shuff = table2array(dataset(dataset.Timepoint=='12','Rshuff_avg')); 
d18shuff = table2array(dataset(dataset.Timepoint=='18','Rshuff_avg')); 
d46shuff = table2array(dataset(dataset.Timepoint=='46','Rshuff_avg'));
alltimesshuff=[BL,d3,d6,d12,d18,d46];

mean_BLshuff = mean(BLshuff, "omitnan"); std_BLshuff = std(BLshuff, "omitnan")/sqrt(length(BLshuff)); % computation of the std error deviation
mean_d3shuff = mean(d3shuff, "omitnan"); std_d3shuff = std(d3shuff, "omitnan")/sqrt(length(d3shuff));
mean_d6shuff = mean(d6shuff, "omitnan"); std_d6shuff = std(d6shuff, "omitnan")/sqrt(length(d6shuff));
mean_d12shuff = mean(d12shuff, "omitnan"); std_d12shuff = std(d12shuff, "omitnan")/sqrt(length(d12shuff));
mean_d18shuff = mean(d18shuff, "omitnan"); std_d18shuff = std(d18shuff, "omitnan")/sqrt(length(d18shuff));
mean_d46shuff = mean(d46shuff, "omitnan"); std_d46shuff = std(d46shuff, "omitnan")/sqrt(length(d46shuff));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];
means_shuff = [mean_BLshuff; mean_d3shuff; mean_d6shuff; mean_d12shuff; mean_d18shuff; mean_d46shuff]; 
stdErrors_shuff = [std_BLshuff; std_d3shuff; std_d6shuff; std_d12shuff; std_d18shuff; std_d46shuff];

% Plot the bar chart with error bars for the standard errors
figure_xorr_D1_allcells=bar(means,'FaceColor',[0 0 1],'FaceAlpha',0.5);
hold on
bar(means_shuff,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.5);

for k = 1:size(means,1)
    scatter(k*(ones(height(alltimes),1)),alltimes(:,k),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','b','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
errorbar(means_shuff, stdErrors_shuff, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('Xcorr - (R-value)'); % Choice of the name of the variable plot HERE
title('Cross Correlation of Spont. Activity of Cells in D1 barrel over time'); % Choice of the name of the plot HERE
ylim([-0.4 1]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for D1 barrel
normality = lillietest(dataset.Rreal_avg) %1=not normal
% lme = fitlme(dataset,'Rreal_avg ~ Timepoint + (1 | MouseID:index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'Rreal_avg ~ Timepoint + (Timepoint | MouseID:index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality aov lme lme2 stdErrors Timepoints
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    BL d3 d6 d12 d18 d46 BLshuff d3shuff d6shuff d12shuff d18shuff d46shuff...
    k means mean_BLshuff mean_d3shuff mean_d6shuff mean_d12shuff mean_d18shuff mean_d46shuff...
    std_BLshuff std_d3shuff std_d6shuff std_d12shuff std_d18shuff std_d46shuff...
    means_shuff stdErrors_shuff alltimes alltimesshuff
%% Fig. S5e: Plot Xcorr of Spont. Activity of cells in C1 barrel over time
dataset = Spont_Data_ByCell_wCorr(Spont_Data_ByCell_wCorr.Barrel_imaged=='C1',:);

dataset=dataset(~isnan(dataset.Rreal_avg),:); %remove cells where xcorr could not be correlated because there was only 1 cell

%Get means for true cross correlation
BL = table2array(dataset(dataset.Timepoint=='0','Rreal_avg')); 
d3 = table2array(dataset(dataset.Timepoint=='3','Rreal_avg')); 
d6 = table2array(dataset(dataset.Timepoint=='6','Rreal_avg'));  
d12 = table2array(dataset(dataset.Timepoint=='12','Rreal_avg')); 
d18 = table2array(dataset(dataset.Timepoint=='18','Rreal_avg')); 
d46 = table2array(dataset(dataset.Timepoint=='46','Rreal_avg'));
alltimes=[BL,d3,d6,d12,d18,d46];

mean_BL = mean(BL, "omitnan"); std_BL = std(BL, "omitnan")/sqrt(length(BL)); % computation of the std error deviation
mean_d3 = mean(d3, "omitnan"); std_d3 = std(d3, "omitnan")/sqrt(length(d3));
mean_d6 = mean(d6, "omitnan"); std_d6 = std(d6, "omitnan")/sqrt(length(d6));
mean_d12 = mean(d12, "omitnan"); std_d12 = std(d12, "omitnan")/sqrt(length(d12));
mean_d18 = mean(d18, "omitnan"); std_d18 = std(d18, "omitnan")/sqrt(length(d18));
mean_d46 = mean(d46, "omitnan"); std_d46 = std(d46, "omitnan")/sqrt(length(d46));

%Get the means for shuffled cross correlation
BLshuff = table2array(dataset(dataset.Timepoint=='0','Rshuff_avg')); 
d3shuff = table2array(dataset(dataset.Timepoint=='3','Rshuff_avg')); 
d6shuff = table2array(dataset(dataset.Timepoint=='6','Rshuff_avg'));  
d12shuff = table2array(dataset(dataset.Timepoint=='12','Rshuff_avg')); 
d18shuff = table2array(dataset(dataset.Timepoint=='18','Rshuff_avg')); 
d46shuff = table2array(dataset(dataset.Timepoint=='46','Rshuff_avg'));
alltimesshuff=[BL,d3,d6,d12,d18,d46];

mean_BLshuff = mean(BLshuff, "omitnan"); std_BLshuff = std(BLshuff, "omitnan")/sqrt(length(BLshuff)); % computation of the std error deviation
mean_d3shuff = mean(d3shuff, "omitnan"); std_d3shuff = std(d3shuff, "omitnan")/sqrt(length(d3shuff));
mean_d6shuff = mean(d6shuff, "omitnan"); std_d6shuff = std(d6shuff, "omitnan")/sqrt(length(d6shuff));
mean_d12shuff = mean(d12shuff, "omitnan"); std_d12shuff = std(d12shuff, "omitnan")/sqrt(length(d12shuff));
mean_d18shuff = mean(d18shuff, "omitnan"); std_d18shuff = std(d18shuff, "omitnan")/sqrt(length(d18shuff));
mean_d46shuff = mean(d46shuff, "omitnan"); std_d46shuff = std(d46shuff, "omitnan")/sqrt(length(d46shuff));

figure
xLabels = {'BL', 'D3', 'D6', 'D12', 'D18', 'D46'}; 
means = [mean_BL; mean_d3; mean_d6; mean_d12; mean_d18; mean_d46]; stdErrors = [std_BL; std_d3; std_d6; std_d12; std_d18; std_d46];
means_shuff = [mean_BLshuff; mean_d3shuff; mean_d6shuff; mean_d12shuff; mean_d18shuff; mean_d46shuff]; 
stdErrors_shuff = [std_BLshuff; std_d3shuff; std_d6shuff; std_d12shuff; std_d18shuff; std_d46shuff];

% Plot the bar chart with error bars for the standard errors
figure_xcorr_C1_allcells=bar(means,'FaceColor',[1 0 1],'FaceAlpha',0.5);
hold on
bar(means_shuff,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.5);

for k = 1:size(means,1)
    scatter(k*(ones(height(alltimes),1)),alltimes(:,k),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','m','XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
errorbar(means_shuff, stdErrors_shuff, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10);
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 6.5])
ylabel('Xcorr - (R-value)'); % Choice of the name of the variable plot HERE
title('Cross Correlation of Spont. Activity of Cells in C1 barrel over time'); % Choice of the name of the plot HERE
ylim([-0.4 1]) % add a limit here if needed to scale different plots together reliability 

%Get statistics for C1 barrel 
normality = lillietest(dataset.Rreal_avg) %1=not normal
% lme = fitlme(dataset,'Rreal_avg ~ Timepoint + (1 | MouseID:index)','FitMethod','REML') %random intercept, cells nested in mice
lme2 = fitlme(dataset,'Rreal_avg ~ Timepoint + (Timepoint | MouseID:index)','FitMethod','REML') %random intercept & slope
% compare(lme,lme2) %adding random slopes does improve the model
aov=anova(lme2)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme2.Coefficients(:,6)))

clear h adj_p adj_ci_cvrg crit_p normality aov lme2 stdErrors
clear dataset xLabels mean_BL mean_d3 mean_d6 mean_d12 mean_d18 mean_d46...
    std_BL std_d3 std_d6 std_d12 std_d18 std_d46...
    BL d3 d6 d12 d18 d46 BLshuff d3shuff d6shuff d12shuff d18shuff d46shuff...
    k means mean_BLshuff mean_d3shuff mean_d6shuff mean_d12shuff mean_d18shuff mean_d46shuff...
    std_BLshuff std_d3shuff std_d6shuff std_d12shuff std_d18shuff std_d46shuff...
    means_shuff stdErrors_shuff alltimes alltimesshuff
%% Fig. S6a: Plot MI of cells in D1 barrel to D1 stim that switch responsivity from BL to D12
dataset = D1_Stm_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12)
NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI_BL_D1(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));
figure
xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI_BL_D1(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI_BL_D1(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in D1 barrel to D1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability 

normality = lillietest(dataset.MI_BL_D1) 
[p,tbl,stats] = kruskalwallis(dataset.MI_BL_D1,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means p stats tbl normality c
%% Fig. S6b: Plot MI of cells in D1 barrel to C1 stim that switch responsivity from BL to D12
dataset = D1_Stm_C1_BL;

NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));

xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_C1_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_C1_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in D1 barrel to C1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability

normality = lillietest(dataset.MI) 
[p,tbl,stats] = kruskalwallis(dataset.MI,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means p stats tbl normality c
%% Fig. S6c: Plot SI of cells in D1 barrel that switch responsivity from BL to D12

dataset = D1_Stm_D1;

NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.SI_MI_BL(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));

xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_SI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_SI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.SI_MI_BL(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.SI_MI_BL(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('SI (C1 vs D1)'); % Choice of the name of the variable plot HERE
title('Selectivity Index of Cells in D1 barrel to D1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-1.1 1.1]) % add a limit here if needed to scale different plots together reliability

normality = lillietest(dataset.SI_MI_BL) 
[p,tbl,stats] = kruskalwallis(dataset.SI_MI_BL,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset NRtoNR NRtoR RtoNR RtoR mean_NRtoNR mean_NRtoR mean_RtoNR...
    mean_RtoR std_NRtoNR std_NRtoR std_RtoNR std_RtoR xLabels means colors...
    c k normality p stats stdErrors SwitchCats tbl
%% Fig. S6d: Plot MI of cells in C1 barrel to D1 stim that switch responsivity from BL to D12

dataset = C1_Stm_D1;

switch_counts=countcats(dataset.Switch_BL_to_D12);

NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI_BL_D1(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI_BL_D1(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));

xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI_BL_D1(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI_BL_D1(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in C1 barrel to D1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability

normality = lillietest(dataset.MI_BL_D1) 
[p,tbl,stats] = kruskalwallis(dataset.MI_BL_D1,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means normality p tbl stats c k switch_counts
%% Fig. S6e: Plot MI of cells in C1 barrel to C1 stim that switch responsivity from BL to D12
dataset = C1_Stm_C1_BL;

NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.MI(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));

xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MI_C1_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MI_C1_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.MI(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.MI(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of Cells in C1 barrel to C1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-5 60]) % add a limit here if needed to scale different plots together reliability

normality = lillietest(dataset.MI) 
[p,tbl,stats] = kruskalwallis(dataset.MI,dataset.Switch_BL_to_D12)
c = multcompare(stats);

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means  normality p tbl stats c k
%% Fig. S6f: Plot SI of cells in C1 barrel that switch responsivity from BL to D12

dataset = C1_Stm_D1;
dataset.SI_MI_BL=C1_BL_D46_all.SI_MI_BL;

NRtoNR = []; NRtoR = []; RtoNR = []; RtoR = [];

for row=1:height(dataset)
    if dataset.Switch_BL_to_D12(row) == "NR to NR" 
        NRtoNR = [NRtoNR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "NR to R" 
        NRtoR = [NRtoR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to NR" 
        RtoNR = [RtoNR dataset.SI_MI_BL(row)]; % change the variable
    elseif dataset.Switch_BL_to_D12(row) == "R to R" 
        RtoR = [RtoR dataset.SI_MI_BL(row)]; % change the variable
    end    
end
clear row

mean_NRtoNR = mean(NRtoNR, "omitnan"); std_NRtoNR = std(NRtoNR, "omitnan")/sqrt(width(NRtoNR)); % computation of the std error deviation
mean_NRtoR = mean(NRtoR, "omitnan"); std_NRtoR = std(NRtoR, "omitnan")/sqrt(width(NRtoR));
mean_RtoR = mean(RtoR, "omitnan"); std_RtoR = std(RtoR, "omitnan")/sqrt(width(RtoR));
mean_RtoNR = mean(RtoNR, "omitnan"); std_RtoNR = std(RtoNR, "omitnan")/sqrt(width(RtoNR));

xLabels = {'NR to NR', 'NR to R', 'R to NR', 'R to R'}; 
means = [mean_NRtoNR; mean_NRtoR; mean_RtoNR; mean_RtoR]; stdErrors = [std_NRtoNR; std_NRtoR; std_RtoNR; std_RtoR];
colors=[0.3 0.3 0.3; 0.7 0.7 0.7; 1 0.7 0.7; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_SI_BL_D12=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_SI_BL_D12.CData(k,:) = colors(k,:);
end
SwitchCats=unique(dataset.Switch_BL_to_D12);
for k = 1:size(SwitchCats,1)
    scatter(k*(ones(length(dataset.SI_MI_BL(dataset.Switch_BL_to_D12==SwitchCats(k))),1)),...
        dataset.SI_MI_BL(dataset.Switch_BL_to_D12==SwitchCats(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlim([0.5 4.5])
xlabel('Switch Category'); 
ylabel('SI (C1 vs D1)'); % Choice of the name of the variable plot HERE
title('Selectivity Index of Cells in C1 barrel to D1 stim at baseline'); % Choice of the name of the plot HERE
ylim([-1.1 1.1]) % add a limit here if needed to scale different plots together reliability

normality = lillietest(dataset.SI_MI_BL) 
[p,tbl,stats] = kruskalwallis(dataset.SI_MI_BL,dataset.Switch_BL_to_D12)
c = multcompare(stats);

C1_Stm_D1 = dataset;

clear dataset xLabels SwitchCats switch_counts NRtoNR NRtoR RtoNR RtoR std_NRtoNR...
    std_NRtoR std_RtoNR std_RtoR stdErrors mean_NRtoNR mean_NRtoR mean_RtoNR mean_RtoR...
    i k colors means normality p tbl stats c k

%% Fig. S6g-h: Perform PCA and GLME analyses to differentiate NR cells that switch responsivity from BL to D12
%First create a datatable with response (switch categories) and predictors
%(Barrel Imaged, MI to C1, MI to D1, SI)
%Get/organize data for C1 barrel
C1_Stm_C1_BL=sortrows(C1_Stm_C1_BL,"Index");
C1_Stm_D1=sortrows(C1_Stm_D1,"Index");
data_C1barrel=C1_Stm_C1_BL(:,{'MouseID','Barrel_imaged','Index','MI','Switch_BL_to_D12','SI_MI_BL'});
data_C1barrel.MI_D1stim=C1_Stm_D1.MI_BL_D1;
data_C1barrel=renamevars(data_C1barrel,['MI'],['MI_C1stim']);
data_C1barrel=movevars(data_C1barrel,'Switch_BL_to_D12','After','SI_MI_BL');
data_C1barrel=movevars(data_C1barrel,'MI_D1stim','Before','SI_MI_BL');
switch_counts_C1=countcats(data_C1barrel.Switch_BL_to_D12);
%Get/organize data for D1 barrel
D1_Stm_D1=sortrows(D1_Stm_D1,"Index");
data_D1barrel=D1_Stm_D1(:,{'MouseID','Barrel_imaged','Index','MI_BL_D1','Switch_BL_to_D12','SI_MI_BL'});
data_D1barrel=renamevars(data_D1barrel,["MI_BL_D1"],["MI_D1stim"]);
data_D1barrel.MI_C1stim=D1_Stm_C1_BL.MI;
data_D1barrel=movevars(data_D1barrel,'MI_C1stim','Before','MI_D1stim');
data_D1barrel=movevars(data_D1barrel,'Switch_BL_to_D12','After','SI_MI_BL');
switch_counts_D1=countcats(data_D1barrel.Switch_BL_to_D12);
%Combine datasetes
dataset=[data_C1barrel;data_D1barrel];
switch_counts_all=countcats(dataset.Switch_BL_to_D12)
%Organize categories so that NR to R is the reference
dataset.Switch_BL_to_D12=reordercats(dataset.Switch_BL_to_D12,...
    {'R to R','R to NR','NR to NR','NR to R'});

%Add in spontaneous data from baseline as a predictor
%Sort rows of spontaneous data to match dataset from above
Spont_BL_Data=Spont_Data_ByCell(Spont_Data_ByCell.Timepoint=='0',:);
Spont_BL_Data=sortrows(Spont_BL_Data,{'Barrel_imaged','MouseID'});
%Add in baseline spontaneous activity metric
dataset.Spont_mean_Z=Spont_BL_Data.mean_Z;

%Make a simplified dataset, keeping only 'NR to NR' vs 'NR to R'
dataset2=dataset(dataset.Switch_BL_to_D12=='NR to NR' | dataset.Switch_BL_to_D12=='NR to R',:);
dataset2.Switch_BL_to_D12=removecats(dataset2.Switch_BL_to_D12,...
    {'R to NR','R to R'});
dataset2.Switch_BL_to_D12=reordercats(dataset2.Switch_BL_to_D12,...
    {'NR to NR','NR to R'});

%Run PCA to determine if cells can be separated by Switch Category
pca_dataset=dataset2;
pca_dataset=sortrows(pca_dataset,'Switch_BL_to_D12');
cats=unique(pca_dataset.Switch_BL_to_D12);
switch_counts_pca=countcats(pca_dataset.Switch_BL_to_D12)
pca_data=[double(pca_dataset.Barrel_imaged),pca_dataset.MI_C1stim,pca_dataset.MI_D1stim,pca_dataset.SI_MI_BL,pca_dataset.Spont_mean_Z];
[coeff,score,latent,tsquared,explained] = pca(pca_data)
explained %First two PCs (PC1>>PC2) explain almost all the variance

s=scatter(score(1:switch_counts_pca(1),1),score(1:switch_counts_pca(1),2),...
    [],[0.3 0.3 0.3],'filled');
hold on
scatter(gca,score((switch_counts_pca(1)+1):end,1),...
    score((switch_counts_pca(1)+1):end,2),[],[0.7 0.7 0.7],'filled');
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
legend(cats)
hold off

%Now fit a glme binomial model to the data
dataset2.Switch_BL_to_D12=double(dataset2.Switch_BL_to_D12);
dataset2.Switch_BL_to_D12(dataset2.Switch_BL_to_D12==2)=0; %NRtoR=0; NRtoNR=1;
glme=fitglme(dataset2,...
    'Switch_BL_to_D12 ~ Barrel_imaged + MI_C1stim + MI_D1stim + SI_MI_BL + Spont_mean_Z + (1|MouseID)',...
    'Distribution','Binomial')
scores_glme=glme.Fitted;

%Perform ROC analysis
[X,Y,T,AUC] = perfcurve(dataset2.Switch_BL_to_D12,scores_glme,...
    '1','NBoot',1000,'XVals',[0:0.01:1]);
AUC %AUC=~0.77 w/ or w/o spont data
figure
e=errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1))
hold on
lx=[0:0.1:1];
ly=lx;
plot(lx,ly)
xlim([-0.05 1.05])
ylim([-0.05 1.05])
xlabel('False positive rate')
ylabel('True positive rate')

clear dataset dataset2 pca_dataset pca_data glme s coeff score latent...
    tsquared explained scores_glme X Y T AUC e lx ly switch_counts_all...
    switch_counts_C1 switch_counts_D1 switch_counts_pca cats data_C1barrel...
    data_D1barrel data_core Spont_BL_Data

%% Fig. 3a: Plot change in responsivity of cells in D1 barrel to D1 stim - BL vs D12 vs D46
Sandata=D1_Stm_D1(:,[12,18,17]);   %need to compute switches @ D12
Sandata.Responsiveness_BL_D1=categorical(Sandata.Responsiveness_BL_D1);
Sandata.Responsiveness_D46_D1=categorical(Sandata.Responsiveness_D46_D1);
Sandata.Responsiveness_BL_D1=renamecats(Sandata.Responsiveness_BL_D1,{'Non-responsive','Responsive'});
Sandata.Responsiveness_D46_D1=renamecats(Sandata.Responsiveness_D46_D1,{'Non-responsive','Responsive'});
Sandata.Switch_BL_to_D12=reordercats(Sandata.Switch_BL_to_D12,{'NR to NR','R to NR', 'NR to R', 'R to R'});
Sandata.Switch_BL_to_D12=renamecats(Sandata.Switch_BL_to_D12,...
    {'NR to NR','R to NR', 'NR to R', 'R to R'},...
    {'_NR to NR','_R to NR', 'NR to R', 'R to R'});
Sandata=sortrows(Sandata,'Switch_BL_to_D12');

options.color_map = [0.3 0.3 0.3;1 0 0;1 0 0;1 0 0;...
    0.3 0.3 0.3;0.3 0.3 0.3;1 0 0;1 0 0;...
    0.3 0.3 0.3;1 0 0;1 0 0;1 0 0;];   
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = false;          % show percentage over the blocks
options.text_color = [1 1 1];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;
plotSankeyFlowChart(Sandata,options);
title('Change in Responsivity of Cells in D1 barrel to D1 - BL vs D12 vs D46');

clear Sandata options
%% Fig. 3b: Plot change in responsivity of cells in C1 barrel to D1 stim - BL vs D12 vs D46
Sandata=C1_Stm_D1(:,[12,18,17]);   
Sandata.Responsiveness_BL_D1=categorical(Sandata.Responsiveness_BL_D1);
Sandata.Responsiveness_D46_D1=categorical(Sandata.Responsiveness_D46_D1);
Sandata.Responsiveness_BL_D1=renamecats(Sandata.Responsiveness_BL_D1,{'Non-responsive','Responsive'});
Sandata.Responsiveness_D46_D1=renamecats(Sandata.Responsiveness_D46_D1,{'Non-responsive','Responsive'});
Sandata.Switch_BL_to_D12=reordercats(Sandata.Switch_BL_to_D12,{'NR to NR','R to NR', 'NR to R', 'R to R'});
Sandata.Switch_BL_to_D12=renamecats(Sandata.Switch_BL_to_D12,...
    {'NR to NR','R to NR', 'NR to R', 'R to R'},...
    {'_NR to NR','_R to NR', 'NR to R', 'R to R'});
Sandata=sortrows(Sandata,'Switch_BL_to_D12');

options.color_map = [0.3 0.3 0.3;1 0 0;1 0 0;1 0 0;...
    0.3 0.3 0.3;0.3 0.3 0.3;1 0 0;1 0 0;...
    0.3 0.3 0.3;1 0 0;1 0 0;1 0 0;];   
options.flow_transparency = 0.2;   % opacity of the flow paths
options.bar_width = 40;            % width of the category blocks
options.show_perc = false;          % show percentage over the blocks
options.text_color = [1 1 1];      % text color for the percentages
options.show_layer_labels = true;  % show layer names under the chart
options.show_cat_labels = true;   % show categories over the blocks.
options.show_legend = false;
plotSankeyFlowChart(Sandata,options);
title('Change in Responsivity of Cells in C1 barrel to D1 - BL vs D12 vs D46');

clear Sandata options

%% FIg. 3c: Plot SI at BL vs D46 for the overall population in D1 barrel using MI

% Choice of the dataset HERE
dataset = D1_BL_D46;

histogram(dataset.SI_MI_BL, 20,"FaceColor", [0.3 0.3 0.3], "FaceAlpha",0.5)
hold on
histogram(dataset.SI_MI_D46, 20 ,"FaceColor", 'b', "FaceAlpha",0.5)
ylim([0 75])
BL_mean = mean(dataset.SI_MI_BL); D46_mean = mean(dataset.SI_MI_D46); yLimits = ylim;
BL_sem = std(dataset.SI_MI_BL)/sqrt(length(dataset.SI_MI_BL));
D46_sem = std(dataset.SI_MI_D46)/sqrt(length(dataset.SI_MI_D46));

line([BL_mean, BL_mean], yLimits, 'Color', [0.3 0.3 0.3], 'LineStyle',':', 'LineWidth', 5); 
line([D46_mean, D46_mean], yLimits, 'Color', 'b', 'LineStyle',':', 'LineWidth', 5);
hold off

xlim([-1.1 1.1])
legend({'BL', 'D46', 'Mean BL', 'Mean D46'})
title('SI in D1 barrel BL v D46')

%Get statistics
data_stats=[D1_BL_D46(:,{'MouseID','Index'});D1_BL_D46(:,{'MouseID','Index'})];
for i=1:size(D1_BL_D46,1)
    data_stats.Timepoint{i}='BL';
    data_stats.Timepoint{i+size(D1_BL_D46,1)}='D46';
end
data_stats.Timepoint=categorical(data_stats.Timepoint);
data_stats.SI=[D1_BL_D46{:,'SI_MI_BL'};D1_BL_D46{:,'SI_MI_D46'}];

normality=lillietest(data_stats.SI) %1, data is not normal
% lme = fitlme(data_stats,'SI ~ Timepoint + (1| MouseID:Index)') %random intercept w/nested cells
lme2 = fitlme(data_stats,'SI ~ Timepoint + (Timepoint| MouseID:Index)') %random slope & intercept
% compare(lme,lme2) %adding random slope marginally improves the model (p=0.07), keeping random slope for consistency w/C1
clear lme lme2 normality data_stats
clear dataset BL_mean D46_mean yLimits BL_sem D46_sem i

%% Fig. 3d: Plot SI at BL vs D46 for the overall population in C1 barrel using MI
dataset = C1_BL_D46_all;

histogram(dataset.SI_MI_BL, 20,"FaceColor", [0.3 0.3 0.3], "FaceAlpha",0.5)
hold on
histogram(dataset.SI_MI_D46, 20 ,"FaceColor", "m", "FaceAlpha",0.5)
ylim([0 40])

BL_mean = mean(dataset.SI_MI_BL); D46_mean = mean(dataset.SI_MI_D46); yLimits = ylim;
BL_sem = std(dataset.SI_MI_BL)/sqrt(length(dataset.SI_MI_BL));
D46_sem = std(dataset.SI_MI_D46)/sqrt(length(dataset.SI_MI_D46));

line([BL_mean, BL_mean], yLimits, 'Color', [0.3 0.3 0.3], 'LineStyle',':', 'LineWidth', 5);
line([D46_mean, D46_mean], yLimits, 'Color', 'm', 'LineStyle',':', 'LineWidth', 5); 
hold off

xlim([-1.1 1.1])
legend({'BL', 'D46', 'Mean BL', 'Mean D46'},'location','northwest')
title('SI in C1 barrel BL v D46')

%Get statistics
data_stats=[C1_BL_D46_all(:,{'MouseID','Index'});C1_BL_D46_all(:,{'MouseID','Index'})];
for i=1:size(C1_BL_D46_all,1)
    data_stats.Timepoint{i}='BL';
    data_stats.Timepoint{i+size(C1_BL_D46_all,1)}='D46';
end
data_stats.Timepoint=categorical(data_stats.Timepoint);
data_stats.SI=[C1_BL_D46_all{:,'SI_MI_BL'};C1_BL_D46_all{:,'SI_MI_D46'}];

normality=lillietest(data_stats.SI) %1, data is not normal
% lme = fitlme(data_stats,'SI ~ Timepoint + (1| MouseID:Index)') %random intercept w/nested cells
lme2 = fitlme(data_stats,'SI ~ Timepoint + (Timepoint| MouseID:Index)') %random slope & intercept
% compare(lme,lme2) %adding random slope does improve the model
clear lme lme2 normality data_stats
clear dataset BL_mean D46_mean yLimits BL_sem D46_sem i
%% Fig. 3e (left panel): Plot the change in % responsive cells to each whisker in the D1 barrel between BL & d46 
dataset_BL = mean_ByMouse(mean_ByMouse.Timepoint=='BL' & mean_ByMouse.Barrel_imaged=='D1',:); 
dataset_d46 = mean_ByMouse(mean_ByMouse.Timepoint=='Day46' & mean_ByMouse.Barrel_imaged=='D1',:); 

delta_C1 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_C1')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_C1'); 
delta_D1 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_D1')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_D1');  
deltas={delta_C1,delta_D1};
mean_C1=mean(delta_C1);
mean_D1=mean(delta_D1);

std_C1 = std(delta_C1,1)/sqrt(size(delta_C1,1));
std_D1 = std(delta_D1,1)/sqrt(size(delta_D1,1));

means = [mean_C1; mean_D1]; stdErrors = [std_C1; std_D1];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1', 'B2', 'E3'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_PrctResp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none'); % parameters to modify for the design (if needed)
for k = 1:2
    scatter(k*(ones(length(deltas{k}),1)),...
        deltas{k},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('Change in % Responsive Cells'); % Choice of the name of the variable plot HERE
title('Change in % Whisker Responsive Cells in D1 barrel after whisker regrowth'); % Choice of the name of the plot HERE
ylim([-80 40]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset = mean_ByMouse(mean_ByMouse.Barrel_imaged=='D1' & (mean_ByMouse.Timepoint=='BL'|mean_ByMouse.Timepoint=='Day46'),:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim, {'Stm_B2', 'Stm_E3'}); 
dataset.Timepoint = removecats(dataset.Timepoint, {'Day3', 'Day6', 'Day12', 'Day18'}); 
normality = lillietest(dataset.NumRespCells) %0, normal distribution
glme = fitglme(dataset,'NumRespCells ~ Whisker_stim + Timepoint + Timepoint*Whisker_stim + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells,'FitMethod','REMPL')
aov = anova(glme)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))

clear delta_C1 delta_D1 colors
clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means deltas i k
clear dataset normality glme aov h crit_p adj_ci_cvrg adj_p dataset_BL dataset_d46
%% Fig. 3e (right panel): Plot the change in % responsive cells to each whisker in the C1 barrel between BL & d46
dataset_BL = mean_ByMouse(mean_ByMouse.Timepoint=='BL' & mean_ByMouse.Barrel_imaged=='C1',:); 
dataset_d46 = mean_ByMouse(mean_ByMouse.Timepoint=='Day46' & mean_ByMouse.Barrel_imaged=='C1',:); 

delta_C1 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_C1')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_C1'); 
delta_D1 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_D1')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_D1'); 
delta_E3 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_E3')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_E3'); 
delta_B2 = dataset_d46.PercentResponsive(dataset_d46.Whisker_stim=='Stm_B2')-dataset_BL.PercentResponsive(dataset_BL.Whisker_stim=='Stm_B2'); 
deltas={delta_C1,delta_D1,delta_B2,delta_E3};
mean_C1=mean(delta_C1);
mean_D1=mean(delta_D1);
mean_E3=mean(delta_E3);
mean_B2=mean(delta_B2);

std_C1 = std(delta_C1,1)/sqrt(size(delta_C1,1));
std_D1 = std(delta_D1,1)/sqrt(size(delta_D1,1));
std_E3 = std(delta_E3,1)/sqrt(size(delta_E3,1));
std_B2 = std(delta_B2,1)/sqrt(size(delta_B2,1));

means = [mean_C1; mean_D1; mean_B2; mean_E3;]; stdErrors = [std_C1; std_D1; std_B2; std_E3;];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1', 'B2', 'E3'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_PrctResp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none'); % parameters to modify for the design (if needed)
for k = 1:size(deltas,2)
    scatter(k*(ones(length(deltas{k}),1)),...
        deltas{k},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('Change in % Responsive Cells'); % Choice of the name of the variable plot HERE
title('Change in % Whisker Responsive Cells in C1 barrel after whisker regrowth'); % Choice of the name of the plot HERE
ylim([-80 40]) % add a limit here if needed to scale different plots together reliability

%Get statistics
dataset = mean_ByMouse(mean_ByMouse.Barrel_imaged=='C1' & (mean_ByMouse.Timepoint=='BL'|mean_ByMouse.Timepoint=='Day46'),:);
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'}); 
dataset.Timepoint = removecats(dataset.Timepoint, {'Day3', 'Day6', 'Day12', 'Day18'}); 
normality = lillietest(dataset.NumRespCells) 
glme = fitglme(dataset,'NumRespCells ~ Whisker_stim + Timepoint + Timepoint*Whisker_stim + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset.TotalCells,'FitMethod','REMPL')
aov = anova(glme)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))

clear delta_B2 delta_C1 delta_D1 delta_E3 colors
clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means deltas i k
clear dataset normality glme aov h crit_p adj_ci_cvrg adj_p dataset_BL dataset_d46

%% Fig. 3f (left panel): Plot the change in MI to each whisker in the D1 barrel between BL & d46

dataset_BL = Data(Data.Timepoint=='BL' & Data.Barrel_imaged=='D1',:); 
dataset_d46 = Data(Data.Timepoint=='Day46' & Data.Barrel_imaged=='D1',:); 
dataset_BL = sortrows(dataset_BL,{'Index','Whisker_stim'}); %Make sure tables are sorted the same
dataset_d46 = sortrows(dataset_d46,{'Index','Whisker_stim'});

delta_C1 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_C1')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_C1'); 
delta_D1 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_D1')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_D1'); 
deltas={delta_C1,delta_D1};
mean_C1=mean(delta_C1);
mean_D1=mean(delta_D1);

std_C1 = std(delta_C1,1)/sqrt(size(delta_C1,1));
std_D1 = std(delta_D1,1)/sqrt(size(delta_D1,1));

means = [mean_C1; mean_D1]; stdErrors = [std_C1; std_D1];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_dMI=bar(means,'FaceColor','flat');
hold on
for k = 1:size(means,1)
    figure_dMI.CData(k,:) = colors(k,:);
end
for k = 1:2
    scatter(k*(ones(length(deltas{k}),1)),...
        deltas{k},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('Change in MI (Z-score)'); % Choice of the name of the variable plot HERE
title('Change in MI in all cells in D1 barrel after whisker regrowth'); % Choice of the name of the plot HERE
ylim([-50 80]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])
clear mean_C1 mean_D1 std_C1 std_D1 xLabels stdErrors means

%Get statistics
dataset = Data(Data.Barrel_imaged=='D1' & (Data.Timepoint=='BL'|Data.Timepoint=='Day46'),:);
dataset.Whisker_stim = removecats(dataset.Whisker_stim, {'Stm_B2', 'Stm_E3'}); 
dataset.Timepoint = removecats(dataset.Timepoint, {'Day3', 'Day6', 'Day12', 'Day18'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + Timepoint + Timepoint*Whisker_stim + (1 | MouseID)',...
   'FitMethod','REML')
aov = anova(lme)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means deltas
clear delta_C1 delta_D1
clear dataset normality lme aov h crit_p adj_ci_cvrg adj_p colors k
%% Fig. 3f (right panel): Plot the change in MI to each whisker in the C1 barrel between BL & d46 

dataset_BL = Data(Data.Timepoint=='BL' & Data.Barrel_imaged=='C1',:); 
dataset_d46 = Data(Data.Timepoint=='Day46' & Data.Barrel_imaged=='C1',:); 
dataset_BL = sortrows(dataset_BL,{'Index','Whisker_stim'}); %Make sure tables are sorted the same
dataset_d46 = sortrows(dataset_d46,{'Index','Whisker_stim'});

delta_C1 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_C1')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_C1'); 
delta_D1 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_D1')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_D1'); 
delta_E3 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_E3')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_E3'); 
delta_B2 = dataset_d46.MI(dataset_d46.Whisker_stim=='Stm_B2')-dataset_BL.MI(dataset_BL.Whisker_stim=='Stm_B2');
deltas={delta_C1,delta_D1,delta_B2,delta_E3};

mean_C1=mean(delta_C1);
mean_D1=mean(delta_D1);
mean_E3=mean(delta_E3);
mean_B2=mean(delta_B2);

std_C1 = std(delta_C1,1)/sqrt(size(delta_C1,1));
std_D1 = std(delta_D1,1)/sqrt(size(delta_D1,1));
std_E3 = std(delta_E3,1)/sqrt(size(delta_E3,1));
std_B2 = std(delta_B2,1)/sqrt(size(delta_B2,1));

means = [mean_C1; mean_D1; mean_B2; mean_E3;]; stdErrors = [std_C1; std_D1; std_B2; std_E3;];

% Plot the bar chart with error bars for the standard errors
xLabels = {'C1', 'D1', 'B2', 'E3'};
colors=[1 0 1; 0 0 1; 0 1 0; 1 0.4 0];

figure_dMI=bar(means,'FaceColor','flat');
hold on
for k = 1:size(means,1)
    figure_dMI.CData(k,:) = colors(k,:);
end

for k = 1:size(deltas,2)
    scatter(k*(ones(length(deltas{k}),1)),...
        deltas{k},100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.5);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulated'); 
ylabel('Change in MI (Z-score)'); % Choice of the name of the variable plot HERE
title('Change in MI in all cells in C1 barrel after whisker regrowth'); % Choice of the name of the plot HERE
ylim([-50 80]) % add a limit here if needed to scale different plots together reliability

%Get statistics
dataset = Data(Data.Barrel_imaged=='C1' & (Data.Timepoint=='BL'|Data.Timepoint=='Day46'),:);
dataset.Whisker_stim = reordercats(dataset.Whisker_stim, {'Stm_C1', 'Stm_D1', 'Stm_B2', 'Stm_E3'}); 
dataset.Timepoint = removecats(dataset.Timepoint, {'Day3', 'Day6', 'Day12', 'Day18'}); 
normality = lillietest(dataset.MI) 
lme = fitlme(dataset,'MI ~ Whisker_stim + Timepoint + Timepoint*Whisker_stim + (1 | MouseID)',...
   'FitMethod','REML')
aov = anova(lme)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear mean_C1 mean_D1 mean_E3 mean_B2 std_C1 std_D1 std_B2 std_E3 xLabels stdErrors means deltas
clear delta_B2 delta_C1 delta_D1 delta_E3
clear dataset normality lme aov h crit_p adj_ci_cvrg adj_p colors k
%% Fig. S7: Calculate and plot tuning curve width in C1 barrel at baseline and D46
%TCW in C1_BL_D46 at BL and at D46 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = C1_BL_D46;
% Set negative MI & MI to 0 for calculation of SI or TCW
for row=1:height(C1_BL_D46)
    for col=5:20
        if dataset{row,col} < 0
            dataset{row,col}= 0;
        end
    end
end
clear row col

for row=1:height(dataset)
    % Computation of SI: (MI_C1 - MI_D1) / (MI_C1 + MI_D1)
    % SI belongs to [-1,1]
    % SI = 0 means no preference for either C1 or D1 (same response)
    % SI = -1 means total selectivity for D1
    % SI = 1 means total selectivity for C1
    dataset.SI_BL(row) = (dataset.MI_BL_C1(row) - dataset.MI_BL_D1(row)) / (dataset.MI_BL_C1(row) + dataset.MI_BL_D1(row));
    dataset.SI_D46(row) = (dataset.MI_D46_C1(row) - dataset.MI_D46_D1(row)) / (dataset.MI_D46_D1(row) + dataset.MI_D46_C1(row));
    dataset.SI_MI_BL(row) = (dataset.MI_BL_C1(row) - dataset.MI_BL_D1(row)) / (dataset.MI_BL_C1(row) + dataset.MI_BL_D1(row));
    dataset.SI_MI_D46(row) = (dataset.MI_D46_C1(row) - dataset.MI_D46_D1(row)) / (dataset.MI_D46_D1(row) + dataset.MI_D46_C1(row));
end

%Set SI=NaN (meaning no response to C1 or D1 whisker) to 0
dataset.SI_BL(isnan(dataset.SI_BL))=0;
dataset.SI_D46(isnan(dataset.SI_D46))=0;
dataset.SI_MI_BL(isnan(dataset.SI_MI_BL))=0;
dataset.SI_MI_D46(isnan(dataset.SI_MI_D46))=0;

%Normalize all MI to the preferred whisker for each neuron
for row=1:height(dataset)
    maximum_BL = max([dataset.MI_BL_C1(row) dataset.MI_BL_D1(row) dataset.MI_BL_E3(row) dataset.MI_BL_B2(row)]);
    maximum_D46 = max([dataset.MI_D46_C1(row) dataset.MI_D46_D1(row) dataset.MI_D46_E3(row) dataset.MI_D46_B2(row)]);
    dataset.normalized_MI_BL_B2(row) = dataset.MI_BL_B2(row)/maximum_BL;
    dataset.normalized_MI_BL_C1(row) = dataset.MI_BL_C1(row)/maximum_BL;   
    dataset.normalized_MI_BL_D1(row) = dataset.MI_BL_D1(row)/maximum_BL;
    dataset.normalized_MI_BL_E3(row) = dataset.MI_BL_E3(row)/maximum_BL;
    dataset.normalized_MI_D46_B2(row) = dataset.MI_D46_B2(row)/maximum_D46;
    dataset.normalized_MI_D46_C1(row) = dataset.MI_D46_C1(row)/maximum_D46;   
    dataset.normalized_MI_D46_D1(row) = dataset.MI_D46_D1(row)/maximum_D46;
    dataset.normalized_MI_D46_E3(row) = dataset.MI_D46_E3(row)/maximum_D46;    
    clear maximum_D46 maximum_BL
end
clear row

% function decay exponential to fit on the data with p(1) a non-important constant and p(2) the decay time
fun = @(p,x) p(1)*exp(-p(2)*x);
figure
for row=1:height(dataset)
    % organization of data + check that we don't have NaN data for all the MIs of a timepoint (which creates an error)
    if (isnan(dataset.normalized_MI_BL_C1(row)) & isnan(dataset.normalized_MI_BL_D1(row)) & isnan(dataset.normalized_MI_BL_B2(row)) & isnan(dataset.normalized_MI_BL_E3(row)))==0
        A = sort([dataset.normalized_MI_BL_C1(row) dataset.normalized_MI_BL_D1(row) dataset.normalized_MI_BL_B2(row) dataset.normalized_MI_BL_E3(row)], 'descend');
    else 
        A = [0 0 0 0]; % if all is NaN: gives the zero value 
    end
    % dataset for fitting
    xdata = [1 2 3 4]; % just 4 indexes for the 4 whiskers
    ydata = A;

    % Initial value to solve the constants (normalized on the hightest value i.e the first point is always (0,1)
    x0 = [0 1];

    % Fit the function to the data
    [p,resnorm,residual] = lsqcurvefit(fun,x0,xdata,ydata);
    
    % Registration of the decay constant and the residuals of the fit in our dataset
    dataset.TCW_BL(row) = p(2); % decay time = Tuning Curve Width. Higher TCW = narrower tuning curve; TCW=0 would be flat
    dataset.resnorm_BL(row) = resnorm; % value of the squared 2-norm of the residual at x: sum((fun(x,xdata)-ydata).^2)
    dataset.residual_BL{row,1}= residual; % value of the residual fun(x,xdata)-ydata at the solution x,    
   
    % plots of function
    normalization_factor = fun(p,1);
%     plot(xdata, fun(p,xdata), '-o') % if we show the fit where the curves do not necessary start at 1
    plot(xdata, fun(p,xdata)/normalization_factor, '-k') % if we want to start all the curves at 1

    hold on  
    clear residual resnorm x0 xdata ydata  p  A normalization_factor 
end 


clear row
xticklabels({'Whisker A', '', 'Whisker B', '', 'Whisker C','', 'Whisker D',})
xlabel('Whisker stimulated')
xlim([0.9 4.1]);
ylim([-0.1 1.3]);
ylabel('Fitted decay exponential')
title(sprintf('C1 barrel TCW, normalized, Baseline')) % Title to modify

hold off
figure
for row=1:height(dataset)
    % organization of data + check that we don't have NaN data for all the MIs of a timepoint (which creates an error)
    if (isnan(dataset.normalized_MI_D46_C1(row)) & isnan(dataset.normalized_MI_D46_D1(row)) & isnan(dataset.normalized_MI_D46_B2(row)) & isnan(dataset.normalized_MI_D46_E3(row)))==0
        A = sort([dataset.normalized_MI_D46_C1(row) dataset.normalized_MI_D46_D1(row) dataset.normalized_MI_D46_B2(row) dataset.normalized_MI_D46_E3(row)], 'descend');
    else 
        A = [0 0 0 0]; % if all is NaN: gives the zero value 
    end
    % dataset for fitting
    xdata = [1 2 3 4]; % just 4 indexes for the 4 whiskers
    ydata = A;

    % Initial value to solve the constants (normalized on the hightest value i.e the first point is always (0,1)
    x0 = [0 1];

    % Fit the function to the data
    [p,resnorm,residual] = lsqcurvefit(fun,x0,xdata,ydata);
    
    % Registration of the decay constant and the residuals of the fit in our dataset
    dataset.TCW_D46(row) = p(2); % decay time = Tuning Curve Width. Higher TCW = narrower tuning curve; TCW=0 would be flat
    dataset.resnorm_D46(row) = resnorm; % value of the squared 2-norm of the residual at x: sum((fun(x,xdata)-ydata).^2)
    dataset.residual_D46{row,1}= residual; % value of the residual fun(x,xdata)-ydata at the solution x,    
   
    % plots of function
    normalization_factor = fun(p,1);
%     plot(xdata, fun(p,xdata), '-o') % if we show the fit where the curves do not necessary start at 1
    plot(xdata, fun(p,xdata)/normalization_factor, '-m') % if we want to start all the curves at 1

    hold on  
    clear residual resnorm x0 xdata ydata  p  A normalization_factor 
end 

clear row
xticklabels({'Whisker A', '', 'Whisker B', '', 'Whisker C','', 'Whisker D',})
xlabel('Whisker stimulated')
xlim([0.9 4.1]);
ylim([-0.1 1.3]);
ylabel('Fitted decay exponential')
title(sprintf('C1 barrel TCW, normalized, D46')) % Title to modify
hold off

figure
histogram(dataset.TCW_BL, 10,"FaceColor", [0.3 0.3 0.3], "FaceAlpha",0.5)
hold on
histogram(dataset.TCW_D46, 10,"FaceColor", "m", "FaceAlpha",0.5)

BL_mean = mean(dataset.TCW_BL); D46_mean = mean(dataset.TCW_D46); yLimits = ylim;
BL_sem = std(dataset.TCW_BL)/sqrt(length(dataset.TCW_BL));
D46_sem = std(dataset.TCW_D46)/sqrt(length(dataset.TCW_D46));

line([BL_mean, BL_mean], yLimits, 'Color', [0.3 0.3 0.3], 'LineStyle',':', 'LineWidth', 5);
line([D46_mean, D46_mean], yLimits, 'Color', 'm', 'LineStyle',':', 'LineWidth', 5); 
hold off

legend({'BL', 'D46', 'Mean BL', 'Mean D46'},'location','northeast')
title('TCW in C1 barrel BL v D46')

data_stats=[dataset(:,{'MouseID','Index'});dataset(:,{'MouseID','Index'})];
for i=1:size(dataset,1)
    data_stats.Timepoint{i}='BL';
    data_stats.Timepoint{i+size(dataset,1)}='D46';
end
data_stats.Timepoint=categorical(data_stats.Timepoint);
data_stats.TCW=[dataset{:,'TCW_BL'};dataset{:,'TCW_D46'}];

normality=lillietest(data_stats.TCW) %1, data is not normal
% lme = fitlme(data_stats,'TCW ~ Timepoint + (1| MouseID:Index)') %random
% intercept with nested cells
lme2 = fitlme(data_stats,'TCW ~ Timepoint + (Timepoint| MouseID:Index)') %random intercept & slope
% compare(lme,lme2) %adding random slope improves the model
clear lme lme2 normality data_stats BL_mean BL_sem D46_mean D46_sem i yLimits
%% Section 4 (Required): Load data to analyze Awake vs. Anesthetized comparisons
%First clear any workspace variables 
clear
%Now load the Awake vs. Anesthetized data table
load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/PVG Anesthetized Awake Comparison/AnesthVsAwake_ByCell.mat')

%Get means for each condition/stim from individual mice
Data=AnesthVsAwake_ByCell;
Data.stim_amp=reordercats(Data.stim_amp,{'Amp2_Spd0.001','Amp2_Spd0.01',...
    'Amp2','Amp4','Amp6','Amp8'});
omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);
T1=Data(:,{'MouseID','condition','stim_amp','whisker_stim'});
[G,mean_ByMouse]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','stim_amp','whisker_stim'});

mean_ByMouse.NumRespCells=splitapply(osum,Data.Responsivity,G);
mean_ByMouse.TotalCells=C.GroupCount;
mean_ByMouse.PercentResponsive=(mean_ByMouse.NumRespCells./mean_ByMouse.TotalCells)*100;

%Then summarize by group to get mean with n=# of mice
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse(:,{'condition','stim_amp','whisker_stim'});
[G,mean_ByGroup]=findgroups(T1);

mean_ByGroup.PercentResponsive=splitapply(omean,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.PercentResponsive_SEM=splitapply(sem,mean_ByMouse.PercentResponsive,G);

% Pool all cells by group for responsive cells and non-responsive cells
Data_Resp=Data(Data.Responsivity==1,:);
T1=Data_Resp(:,{'condition','stim_amp','whisker_stim'});
[G,mean_Resp]=findgroups(T1);
C=groupcounts(T1,{'condition','stim_amp','whisker_stim'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_Resp.TotalCells=C.GroupCount;
mean_Resp.AUC=splitapply(omean,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.AUC_SEM=splitapply(sem,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.MI=splitapply(omean,Data_Resp.MI,G);
mean_Resp.MI_SEM=splitapply(sem,Data_Resp.MI,G);
mean_Resp.Amp=splitapply(omean,Data_Resp.PeaksMax_mean,G);
mean_Resp.Amp_SEM=splitapply(sem,Data_Resp.PeaksMax_mean,G);
mean_Resp.Latency=splitapply(omean,Data_Resp.PeaksLatency_mean,G);
mean_Resp.Latency_SEM=splitapply(sem,Data_Resp.PeaksLatency_mean,G);
mean_Resp.PrctStimActive=splitapply(omean,Data_Resp.PrctStimActive,G);
mean_Resp.PrctStimActive_SEM=splitapply(sem,Data_Resp.PrctStimActive,G);
mean_Resp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp.allstimsROIs_mean),G);
mean_Resp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp.allstimsROIs_mean),G);

Data_NResp=Data(Data.Responsivity==0,:);
T1=Data_NResp(:,{'condition','stim_amp','whisker_stim'});
[G,mean_NResp]=findgroups(T1);
C=groupcounts(T1,{'condition','stim_amp','whisker_stim'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_NResp.TotalCells=C.GroupCount;
mean_NResp.AUC=splitapply(omean,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.AUC_SEM=splitapply(sem,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.MI=splitapply(omean,Data_NResp.MI,G);
mean_NResp.MI_SEM=splitapply(sem,Data_NResp.MI,G);
mean_NResp.Amp=splitapply(omean,Data_NResp.PeaksMax_mean,G);
mean_NResp.Amp_SEM=splitapply(sem,Data_NResp.PeaksMax_mean,G);
mean_NResp.Latency=splitapply(omean,Data_NResp.PeaksLatency_mean,G);
mean_NResp.Latency_SEM=splitapply(sem,Data_NResp.PeaksLatency_mean,G);
mean_NResp.PrctStimActive=splitapply(omean,Data_NResp.PrctStimActive,G);
mean_NResp.PrctStimActive_SEM=splitapply(sem,Data_NResp.PrctStimActive,G);
mean_NResp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp.allstimsROIs_mean),G);
mean_NResp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp.allstimsROIs_mean),G);

clear C G omean omean2 osum sem sem2 T1
%% Fig. S1e: Plot the % responsive cells in Awake & Anesthetized conditions by whisker stimulus type
dataset = mean_ByGroup(mean_ByGroup.whisker_stim=='D1',:); 
dataset = sortrows(dataset,'stim_amp');
means= [table2array(dataset(dataset.condition=='Anesth','PercentResponsive')),...
    table2array(dataset(dataset.condition=='Awake','PercentResponsive'))];
stdErrors = [table2array(dataset(dataset.condition=='Anesth','PercentResponsive_SEM')),...
    table2array(dataset(dataset.condition=='Awake','PercentResponsive_SEM'))];
stimtypes=unique(dataset.stim_amp);
conditions=unique(dataset.condition);
% xLabels = {'Amp2', 'Amp4', 'Amp6', 'Amp8','RT 0.01','RT 0.001'};
xLabels=stimtypes;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on
errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none');
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none');

for j= 1:2 %condition
    for k = 1:size(stimtypes,1)
        plotdata=mean_ByMouse.PercentResponsive(mean_ByMouse.condition==conditions(j) & mean_ByMouse.stim_amp==stimtypes(k) & mean_ByMouse.whisker_stim=='D1');
        scatter(xdata(k,j)*(ones(length(plotdata),1)),...
            plotdata,50,'Marker','o',...
            'MarkerEdgeColor','k','MarkerFaceColor',colors(j,:),'XJitter','randn','XJitterWidth',0.2);
    end
end
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Whisker Stimulation'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('% Whisker Responsive Cells by Stimulation Type'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability
legend({'Anesth','Awake'},'Location','NorthWest')

%Get statistics
dataset2=mean_ByMouse(mean_ByMouse.whisker_stim=='D1',:);
normality = lillietest(dataset2.PercentResponsive) %1, not normal
glme = fitglme(dataset2,'NumRespCells ~ stim_amp + condition + condition*stim_amp + (1 | MouseID)',...
   'Distribution','Binomial','BinomialSize', dataset2.TotalCells,'FitMethod','REMPL')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(glme.Coefficients(:,6)))
aov = anova(glme)

clear means StdErrors plotdata xLabels stdErrors means dataset dataset2...
    glme aov normality h j k xdata stimtypes crit_p adj_ci_cvrg adj_p colors conditions

%% Fig. S1f: Plot the mean-evoked traces from responders to the D1 whisker in the D1 barrel, Awake vs. Anesth
dataset = AnesthVsAwake_ByCell(AnesthVsAwake_ByCell.stim_amp=='Amp2_Spd0.001' & ...
    AnesthVsAwake_ByCell.whisker_stim=='D1' & AnesthVsAwake_ByCell.Responsivity==1,:);

Anesth = []; Awake = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "Anesth" 
        Anesth = [Anesth dataset.allstimsROIs_mean(row)]; % change the variable
    elseif dataset.condition(row) == "Awake" 
        Awake = [Awake dataset.allstimsROIs_mean(row)]; % change the variable
    end    
end
clear row

mean_anesth = mean(cell2mat(Anesth'), "omitnan"); std_anesth = std(cell2mat(Anesth'), "omitnan")/sqrt(width(Anesth)); % computation of the std error deviation
mean_awake = mean(cell2mat(Awake'), "omitnan"); std_awake = std(cell2mat(Awake'), "omitnan")/sqrt(width(Awake));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_anesth', mean_awake']; stdErrors = [std_anesth', std_awake'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-k','-r'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_anesth,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('Mean whisker-evoked response of responsive ROIs in C1 barrel at baseline'); % Choice of the name of the plot HERE
ylim([0 20]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend({'','','Anesth','','Awake'},'Location','NorthWest')


clear dataset mean_anesth mean_awake means std_anesth std_awake stdErrors xLabels...
    Anesth Awake greybar i style x

%% Fig. S1g: Plot the MI from Responders in the D1 barrel to the D1 whisker at baseline, Awake vs. Anesth

% Choice of the dataset HERE
dataset = AnesthVsAwake_ByCell(AnesthVsAwake_ByCell.stim_amp=='Amp2_Spd0.001' & ...
    AnesthVsAwake_ByCell.whisker_stim=='D1' & AnesthVsAwake_ByCell.Responsivity==1,:);

Anesth = []; Awake = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "Anesth" 
        Anesth = [Anesth dataset.MI(row)]; % change the variable
    elseif dataset.condition(row) == "Awake" 
        Awake = [Awake dataset.MI(row)]; % change the variable
    end    
end
clear row

mean_anesth = mean(Anesth, "omitnan"); std_anesth = std(Anesth, "omitnan")/sqrt(width(Anesth)); % computation of the std error deviation
mean_awake = mean(Awake, "omitnan"); std_awake = std(Awake, "omitnan")/sqrt(width(Awake));

xLabels = {'Anesth', 'Awake'}; means = [mean_anesth; mean_awake]; stdErrors = [std_anesth; std_awake];
colors=[0.3 0.3 0.3; 1 0 0];
% Plot the bar chart with error bars for the standard errors
figure_MIresp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_MIresp.CData(k,:) = colors(k,:);
end

conditions=unique(dataset.condition);
for k = 1:2
    scatter(k*(ones(length(dataset.MI(dataset.condition==conditions(k))),1)),...
        dataset.MI(dataset.condition==conditions(k)),100,'Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:),'XJitter','randn','XJitterWidth',0.6);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2,'CapSize',10); % parameters to modify for the design (if needed)
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Condition'); 
ylabel('MI - (Z-score)'); % Choice of the name of the variable plot HERE
title('Modulation Index of the responsive ROIs to Amp2-Spd0.001 stimulation'); % Choice of the name of the plot HERE
ylim([-2 20]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
normality = lillietest(dataset.MI) %1, not normal distribution
lme = fitlme(dataset,'MI ~ condition + (1 | MouseID)','FitMethod','REML')
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(double(lme.Coefficients(:,6)))

clear anesth awake mean_anesth mean_awake dataset xLabels stdErrors means...
    adj_ci_cvrg adj_p Anesth Awake colors crit_p h k lme normality std_anesth...
    std_awake conditions
%% Section 5 (Required): Load the Data for PV DREADD Experiments - Analysis of PV and Pyr Cells
clear %clear other workspace variables

load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/Acute_DREADDs_Suite2P/Dose Comparison/Acute_DREADDs_DataByCell_DoseComparison.mat')

% Use the 5mg/kg C21 dose data for primary analyses
Data=Data5; %Can change to Data1 or Data3 to analyze 1 or 3 mg/kg dose data

%% Get means for each condition/group from individual mice and calculate activity metrics using dF/F
omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);
T1=Data(:,{'MouseID','condition','group','type'});
[G,mean_ByMouse]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','group','type'});

mean_ByMouse.NumRespCells=splitapply(osum,Data.Responsivity,G);
mean_ByMouse.TotalCells=C.GroupCount;
mean_ByMouse.PercentResponsive=(mean_ByMouse.NumRespCells./mean_ByMouse.TotalCells)*100;

%Then summarize by group to get mean with n=# of mice
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse(:,{'condition','group','type'});
[G,mean_ByGroup]=findgroups(T1);

mean_ByGroup.PercentResponsive=splitapply(omean,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.PercentResponsive_SEM=splitapply(sem,mean_ByMouse.PercentResponsive,G);

clear C G T1

% Get mean evoked dF_F traces
allstimsallROIs_dF=cell(size(Data,1),1);
framerate=10.0143333333333;
for k=1:size(Data,1)
    stimulus_trace=Data.stim_epochs{k};
    signal_traces=Data.dF_F{k};
    
    idx=find(stimulus_trace);
    stimstarts=idx(1,1);
    stimends=[];
    for i=2:size(idx,1)
        if idx(i)>(idx(i-1)+1)
            stimstarts=[stimstarts,idx(i)];
            stimends=[stimends,idx(i-1)];
        else
        end
    end
    stimends=[stimends,idx(end)];
    stimstarts=stimstarts';
    stimends=stimends';

    numberofstims=size(stimstarts,1);
    timeafterstim = 3.1; %Given an inter-stim interval of 3.5-6 seconds 
    framesafterstim = round(timeafterstim*framerate);
    timebeforestim = 1; %To be used as a pre-stimulus baseline
    framesbeforestim = round(timebeforestim*framerate);
    chopsize = framesbeforestim+framesafterstim+1;
    peakcut=framesbeforestim+round(2*framerate); %Set to baseline plus 2 seconds

    ROIcount = size(signal_traces,1); 

    %creates a 3D matrix of all stimulus-evoked epochs. 3 parameters: ROIs, stimulation (from
    %1 to numberofstims), frames (deltaF data)
    allstimsallROIs = zeros(chopsize,numberofstims);

    %Now indexing into this 3D matrix. Frames x Stims x ROI
    for j = 1:numberofstims %establishes the matrix columns, i.e. the whisker stimulations
        allstimsallROIs(:,j) = (signal_traces((stimstarts(j)-framesbeforestim):(stimstarts(j)+framesafterstim)));
    end  
    allstimsallROIs_dF{k}=allstimsallROIs;
end

for i=1:size(allstimsallROIs_dF,1)
    allstimsallROIs_dF_mean{i}=mean(allstimsallROIs_dF{i},2)';
end
allstimsallROIs_dF_mean=allstimsallROIs_dF_mean';

%Add to the table
Data.allstimsallROIs_dF=allstimsallROIs_dF;
Data.allstimsallROIs_dF_mean=allstimsallROIs_dF_mean;

% Calculate mean and std of baseline F
for ii=1:size(Data,1)
    Fraw=double(Data.F{ii});
    Fneu=double(Data.Fneu{ii});
    F=(Fraw-0.7*Fneu)';
    Fsize=size(F);
    baselineframes=round(10*framerate);

    %-----------Calculate Baseline F------------------------
    basecalc=zeros(1,Fsize(1)-baselineframes+1);
    for i=1:Fsize(1)-baselineframes+1 %find 10 seconds of baseline. RECALC based on framerate
        basecalc(i)=std(F(i:i+baselineframes-1)); %checks STD for all values
    end
    [~,mindex]=min(basecalc); %calculate the position of lowest deviation for baseline
    Data.Fbase_avg(ii)=mean(F(mindex:mindex+baselineframes-1));
    Data.Fbase_std(ii)=basecalc(mindex);
    %----------------End Calculate Baseline F--------------------------------
    clear Fraw Fneu F Fsize baselineframes basecalc mindex
end

Data.Fbase_noise=abs(Data.Fbase_std./Data.Fbase_avg);

% %Calculate MI using dF - can adjust early/late time windows here if needed
% for i=1:height(Data)
%     mean_trace=Data.allstimsallROIs_dF_mean{i};
%     Data.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
%     Data.MIearly_dF(i)=mean(mean_trace(11:15))-mean(mean_trace(1:10));
%     Data.MIlate_dF(i)=mean(mean_trace(26:30))-mean(mean_trace(1:10));
% end

clear allstimsallROIs_dF allstimsallROIs_dF_mean stimulus_trace signal_traces...
    stimends stimstarts timeafterstim timebeforestim peakcut i j k framesbeforestim...
    framesafterstim chopsize idx ROIcount

% Pool all cells by group for responsive cells and non-responsive cells
Data_Resp=Data(Data.Responsivity==1,:);
T1=Data_Resp(:,{'condition','group','type'});
[G,mean_Resp]=findgroups(T1);
C=groupcounts(T1,{'condition','group','type'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,1,'omitnan')/sqrt(size(x,1)));

mean_Resp.TotalCells=C.GroupCount;
mean_Resp.AUC=splitapply(omean,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.AUC_SEM=splitapply(sem,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.MI=splitapply(omean,Data_Resp.MI,G);
mean_Resp.MI_SEM=splitapply(sem,Data_Resp.MI,G);
mean_Resp.Amp=splitapply(omean,Data_Resp.PeaksMax_mean,G);
mean_Resp.Amp_SEM=splitapply(sem,Data_Resp.PeaksMax_mean,G);
mean_Resp.Latency=splitapply(omean,Data_Resp.PeaksLatency_mean,G);
mean_Resp.Latency_SEM=splitapply(sem,Data_Resp.PeaksLatency_mean,G);
mean_Resp.PrctStimActive=splitapply(omean,Data_Resp.PrctStimActive,G);
mean_Resp.PrctStimActive_SEM=splitapply(sem,Data_Resp.PrctStimActive,G);
mean_Resp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp.allstimsROIs_mean),G);
mean_Resp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp.allstimsROIs_mean),G);

Data_NResp=Data(Data.Responsivity==0,:);
T1=Data_NResp(:,{'condition','group','type'});
[G,mean_NResp]=findgroups(T1);
C=groupcounts(T1,{'condition','group','type'});

mean_NResp.TotalCells=C.GroupCount;
mean_NResp.AUC=splitapply(omean,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.AUC_SEM=splitapply(sem,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.MI=splitapply(omean,Data_NResp.MI,G);
mean_NResp.MI_SEM=splitapply(sem,Data_NResp.MI,G);
mean_NResp.Amp=splitapply(omean,Data_NResp.PeaksMax_mean,G);
mean_NResp.Amp_SEM=splitapply(sem,Data_NResp.PeaksMax_mean,G);
mean_NResp.Latency=splitapply(omean,Data_NResp.PeaksLatency_mean,G);
mean_NResp.Latency_SEM=splitapply(sem,Data_NResp.PeaksLatency_mean,G);
mean_NResp.PrctStimActive=splitapply(omean,Data_NResp.PrctStimActive,G);
mean_NResp.PrctStimActive_SEM=splitapply(sem,Data_NResp.PrctStimActive,G);
mean_NResp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp.allstimsROIs_mean),G);
mean_NResp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp.allstimsROIs_mean),G);

clear C G omean omean2 osum sem sem2

% Calculate MI for early (1-2s) and late time windows (2-3s) after stim onset
for i=1:height(Data)
    mean_trace=Data.allstimsallROIs_dF_mean{i};
    Data.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
    Data.MIearly_dF(i)=mean(mean_trace(11:20))-mean(mean_trace(1:10));
    Data.MIlate_dF(i)=mean(mean_trace(21:30))-mean(mean_trace(1:10));
end

clear x i mean_trace coefficients allstimsallROIs framerate ii numberofstims T1

% Get means for each condition/group from individual mice in 1 mg/kg dose
% Calculate MI for early (1-2s) and late time windows (2-3s) after stim onset
for i=1:height(Data1)
    mean_trace=Data1.allstimsallROIs_dF_mean{i};
    Data1.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
    Data1.MIearly_dF(i)=mean(mean_trace(11:20))-mean(mean_trace(1:10));
    Data1.MIlate_dF(i)=mean(mean_trace(21:30))-mean(mean_trace(1:10));
end

omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);
T1=Data1(:,{'MouseID','condition','group','type'});
[G,mean_ByMouse1]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','group','type'});

mean_ByMouse1.NumRespCells=splitapply(osum,Data1.Responsivity,G);
mean_ByMouse1.TotalCells=C.GroupCount;
mean_ByMouse1.PercentResponsive=(mean_ByMouse1.NumRespCells./mean_ByMouse1.TotalCells)*100;
mean_ByMouse1.MI_dF=splitapply(omean,Data1.MI_dF,G);
mean_ByMouse1.MIearly_dF=splitapply(omean,Data1.MIearly_dF,G);
mean_ByMouse1.MIlate_dF=splitapply(omean,Data1.MIlate_dF,G);

%Then summarize by group to get mean with n=# of mice
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse1(:,{'condition','group','type'});
[G,mean_ByGroup1]=findgroups(T1);

mean_ByGroup1.PercentResponsive=splitapply(omean,mean_ByMouse1.PercentResponsive,G);
mean_ByGroup1.PercentResponsive_SEM=splitapply(sem,mean_ByMouse1.PercentResponsive,G);
mean_ByGroup1.MI_dF=splitapply(omean,mean_ByMouse1.MI_dF,G);
mean_ByGroup1.MI_dF_SEM=splitapply(sem,mean_ByMouse1.MI_dF,G);
mean_ByGroup1.MIearly_dF=splitapply(omean,mean_ByMouse1.MIearly_dF,G);
mean_ByGroup1.MIearly_dF_SEM=splitapply(sem,mean_ByMouse1.MIearly_dF,G);
mean_ByGroup1.MIlate_dF=splitapply(omean,mean_ByMouse1.MIlate_dF,G);
mean_ByGroup1.MIlate_dF_SEM=splitapply(sem,mean_ByMouse1.MIlate_dF,G);

clear C G T1

% Get means for each condition/group from individual mice in 3 mg/kg dose
% Calculate MI for early (1-2s) and late time windows (2-3s) after stim onset
for i=1:height(Data3)
    mean_trace=Data3.allstimsallROIs_dF_mean{i};
    Data3.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
    Data3.MIearly_dF(i)=mean(mean_trace(11:20))-mean(mean_trace(1:10));
    Data3.MIlate_dF(i)=mean(mean_trace(21:30))-mean(mean_trace(1:10));
end
T1=Data3(:,{'MouseID','condition','group','type'});
[G,mean_ByMouse3]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','group','type'});

mean_ByMouse3.NumRespCells=splitapply(osum,Data3.Responsivity,G);
mean_ByMouse3.TotalCells=C.GroupCount;
mean_ByMouse3.PercentResponsive=(mean_ByMouse3.NumRespCells./mean_ByMouse3.TotalCells)*100;
mean_ByMouse3.MI_dF=splitapply(omean,Data3.MI_dF,G);
mean_ByMouse3.MIearly_dF=splitapply(omean,Data3.MIearly_dF,G);
mean_ByMouse3.MIlate_dF=splitapply(omean,Data3.MIlate_dF,G);

%Then summarize by group to get mean with n=# of mice
T1=mean_ByMouse3(:,{'condition','group','type'});
[G,mean_ByGroup3]=findgroups(T1);

mean_ByGroup3.PercentResponsive=splitapply(omean,mean_ByMouse3.PercentResponsive,G);
mean_ByGroup3.PercentResponsive_SEM=splitapply(sem,mean_ByMouse3.PercentResponsive,G);
mean_ByGroup3.MI_dF=splitapply(omean,mean_ByMouse3.MI_dF,G);
mean_ByGroup3.MI_dF_SEM=splitapply(sem,mean_ByMouse3.MI_dF,G);
mean_ByGroup3.MIearly_dF=splitapply(omean,mean_ByMouse3.MIearly_dF,G);
mean_ByGroup3.MIearly_dF_SEM=splitapply(sem,mean_ByMouse3.MIearly_dF,G);
mean_ByGroup3.MIlate_dF=splitapply(omean,mean_ByMouse3.MIlate_dF,G);
mean_ByGroup3.MIlate_dF_SEM=splitapply(sem,mean_ByMouse3.MIlate_dF,G);

clear C G T1

% Get means for each condition/group from individual mice in 5 mg/kg dose
% Calculate MI for early (1-2s) and late time windows (2-3s) after stim onset
for i=1:height(Data5)
    mean_trace=Data5.allstimsallROIs_dF_mean{i};
    Data5.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
    Data5.MIearly_dF(i)=mean(mean_trace(11:20))-mean(mean_trace(1:10));
    Data5.MIlate_dF(i)=mean(mean_trace(21:30))-mean(mean_trace(1:10));
end
T1=Data5(:,{'MouseID','condition','group','type'});
[G,mean_ByMouse5]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','group','type'});

mean_ByMouse5.NumRespCells=splitapply(osum,Data5.Responsivity,G);
mean_ByMouse5.TotalCells=C.GroupCount;
mean_ByMouse5.PercentResponsive=(mean_ByMouse5.NumRespCells./mean_ByMouse5.TotalCells)*100;
mean_ByMouse5.MI_dF=splitapply(omean,Data5.MI_dF,G);
mean_ByMouse5.MIearly_dF=splitapply(omean,Data5.MIearly_dF,G);
mean_ByMouse5.MIlate_dF=splitapply(omean,Data5.MIlate_dF,G);

%Then summarize by group to get mean with n=# of mice
T1=mean_ByMouse5(:,{'condition','group','type'});
[G,mean_ByGroup5]=findgroups(T1);

mean_ByGroup5.PercentResponsive=splitapply(omean,mean_ByMouse5.PercentResponsive,G);
mean_ByGroup5.PercentResponsive_SEM=splitapply(sem,mean_ByMouse5.PercentResponsive,G);
mean_ByGroup5.MI_dF=splitapply(omean,mean_ByMouse5.MI_dF,G);
mean_ByGroup5.MI_dF_SEM=splitapply(sem,mean_ByMouse5.MI_dF,G);
mean_ByGroup5.MIearly_dF=splitapply(omean,mean_ByMouse5.MIearly_dF,G);
mean_ByGroup5.MIearly_dF_SEM=splitapply(sem,mean_ByMouse5.MIearly_dF,G);
mean_ByGroup5.MIlate_dF=splitapply(omean,mean_ByMouse5.MIlate_dF,G);
mean_ByGroup5.MIlate_dF_SEM=splitapply(sem,mean_ByMouse5.MIlate_dF,G);

clear C G T1

% Get change in percent responsive and MI_dF pre/post C21
T1=mean_ByMouse1(:,{'MouseID','group','type'});
[G,delta_ByMouse1]=findgroups(T1);
Dprct1_BL=mean_ByMouse1.PercentResponsive(mean_ByMouse1.condition=='BL');
Dprct1_C21=mean_ByMouse1.PercentResponsive(mean_ByMouse1.condition=='C21');
delta_ByMouse1.Dprct=Dprct1_C21-Dprct1_BL;
DMI_dF1_BL=mean_ByMouse1.MI_dF(mean_ByMouse1.condition=='BL');
DMI_dF1_C21=mean_ByMouse1.MI_dF(mean_ByMouse1.condition=='C21');
delta_ByMouse1.DMI_dF=DMI_dF1_C21-DMI_dF1_BL;
DMIearly_dF1_BL=mean_ByMouse1.MIearly_dF(mean_ByMouse1.condition=='BL');
DMIearly_dF1_C21=mean_ByMouse1.MIearly_dF(mean_ByMouse1.condition=='C21');
delta_ByMouse1.DMIearly_dF=DMIearly_dF1_C21-DMIearly_dF1_BL;
DMIlate_dF1_BL=mean_ByMouse1.MIlate_dF(mean_ByMouse1.condition=='BL');
DMIlate_dF1_C21=mean_ByMouse1.MIlate_dF(mean_ByMouse1.condition=='C21');
delta_ByMouse1.DMIlate_dF=DMIlate_dF1_C21-DMIlate_dF1_BL;
T2=delta_ByMouse1(:,{'group','type'});
[G2,delta_ByGroup1]=findgroups(T2);
delta_ByGroup1.Dprct=splitapply(omean,delta_ByMouse1.Dprct,G2);
delta_ByGroup1.Dprct_SEM=splitapply(sem,delta_ByMouse1.Dprct,G2);
delta_ByGroup1.DMI_dF=splitapply(omean,delta_ByMouse1.DMI_dF,G2);
delta_ByGroup1.DMI_dF_SEM=splitapply(sem,delta_ByMouse1.DMI_dF,G2);
delta_ByGroup1.DMIearly_dF=splitapply(omean,delta_ByMouse1.DMIearly_dF,G2);
delta_ByGroup1.DMIearly_dF_SEM=splitapply(sem,delta_ByMouse1.DMIearly_dF,G2);
delta_ByGroup1.DMIlate_dF=splitapply(omean,delta_ByMouse1.DMIlate_dF,G2);
delta_ByGroup1.DMIlate_dF_SEM=splitapply(sem,delta_ByMouse1.DMIlate_dF,G2);
clear T1 G T2 G2

T1=mean_ByMouse3(:,{'MouseID','group','type'});
[G,delta_ByMouse3]=findgroups(T1);
Dprct3_BL=mean_ByMouse3.PercentResponsive(mean_ByMouse3.condition=='BL');
Dprct3_C21=mean_ByMouse3.PercentResponsive(mean_ByMouse3.condition=='C21');
delta_ByMouse3.Dprct=Dprct3_C21-Dprct3_BL;
DMI_dF3_BL=mean_ByMouse3.MI_dF(mean_ByMouse3.condition=='BL');
DMI_dF3_C21=mean_ByMouse3.MI_dF(mean_ByMouse3.condition=='C21');
delta_ByMouse3.DMI_dF=DMI_dF3_C21-DMI_dF3_BL;
DMIearly_dF3_BL=mean_ByMouse3.MIearly_dF(mean_ByMouse3.condition=='BL');
DMIearly_dF3_C21=mean_ByMouse3.MIearly_dF(mean_ByMouse3.condition=='C21');
delta_ByMouse3.DMIearly_dF=DMIearly_dF3_C21-DMIearly_dF3_BL;
DMIlate_dF3_BL=mean_ByMouse3.MIlate_dF(mean_ByMouse3.condition=='BL');
DMIlate_dF3_C21=mean_ByMouse3.MIlate_dF(mean_ByMouse3.condition=='C21');
delta_ByMouse3.DMIlate_dF=DMIlate_dF3_C21-DMIlate_dF3_BL;
T2=delta_ByMouse3(:,{'group','type'});
[G2,delta_ByGroup3]=findgroups(T2);
delta_ByGroup3.Dprct=splitapply(omean,delta_ByMouse3.Dprct,G2);
delta_ByGroup3.Dprct_SEM=splitapply(sem,delta_ByMouse3.Dprct,G2);
delta_ByGroup3.DMI_dF=splitapply(omean,delta_ByMouse3.DMI_dF,G2);
delta_ByGroup3.DMI_dF_SEM=splitapply(sem,delta_ByMouse3.DMI_dF,G2);
delta_ByGroup3.DMIearly_dF=splitapply(omean,delta_ByMouse3.DMIearly_dF,G2);
delta_ByGroup3.DMIearly_dF_SEM=splitapply(sem,delta_ByMouse3.DMIearly_dF,G2);
delta_ByGroup3.DMIlate_dF=splitapply(omean,delta_ByMouse3.DMIlate_dF,G2);
delta_ByGroup3.DMIlate_dF_SEM=splitapply(sem,delta_ByMouse3.DMIlate_dF,G2);
clear T1 G T2 G2

T1=mean_ByMouse5(:,{'MouseID','group','type'});
[G,delta_ByMouse5]=findgroups(T1);
Dprct5_BL=mean_ByMouse5.PercentResponsive(mean_ByMouse5.condition=='BL');
Dprct5_C21=mean_ByMouse5.PercentResponsive(mean_ByMouse5.condition=='C21');
delta_ByMouse5.Dprct=Dprct5_C21-Dprct5_BL;
DMI_dF5_BL=mean_ByMouse5.MI_dF(mean_ByMouse5.condition=='BL');
DMI_dF5_C21=mean_ByMouse5.MI_dF(mean_ByMouse5.condition=='C21');
delta_ByMouse5.DMI_dF=DMI_dF5_C21-DMI_dF5_BL;
DMIearly_dF5_BL=mean_ByMouse5.MIearly_dF(mean_ByMouse5.condition=='BL');
DMIearly_dF5_C21=mean_ByMouse5.MIearly_dF(mean_ByMouse5.condition=='C21');
delta_ByMouse5.DMIearly_dF=DMIearly_dF5_C21-DMIearly_dF5_BL;
DMIlate_dF5_BL=mean_ByMouse5.MIlate_dF(mean_ByMouse5.condition=='BL');
DMIlate_dF5_C21=mean_ByMouse5.MIlate_dF(mean_ByMouse5.condition=='C21');
delta_ByMouse5.DMIlate_dF=DMIlate_dF5_C21-DMIlate_dF5_BL;
T2=delta_ByMouse5(:,{'group','type'});
[G2,delta_ByGroup5]=findgroups(T2);
delta_ByGroup5.Dprct=splitapply(omean,delta_ByMouse5.Dprct,G2);
delta_ByGroup5.Dprct_SEM=splitapply(sem,delta_ByMouse5.Dprct,G2);
delta_ByGroup5.DMI_dF=splitapply(omean,delta_ByMouse5.DMI_dF,G2);
delta_ByGroup5.DMI_dF_SEM=splitapply(sem,delta_ByMouse5.DMI_dF,G2);
delta_ByGroup5.DMIearly_dF=splitapply(omean,delta_ByMouse5.DMIearly_dF,G2);
delta_ByGroup5.DMIearly_dF_SEM=splitapply(sem,delta_ByMouse5.DMIearly_dF,G2);
delta_ByGroup5.DMIlate_dF=splitapply(omean,delta_ByMouse5.DMIlate_dF,G2);
delta_ByGroup5.DMIlate_dF_SEM=splitapply(sem,delta_ByMouse5.DMIlate_dF,G2);
clear T1 G T2 G2 DMI_dF1_BL DMI_dF1_C21 DMI_dF3_BL DMI_dF3_C21 DMI_dF5_BL DMI_dF5_C21...
    DMIearly_dF1_BL DMIearly_dF1_C21 DMIearly_dF3_BL DMIearly_dF3_C21 DMIearly_dF5_BL DMIearly_dF5_C21...
    DMIlate_dF1_BL DMIlate_dF1_C21 DMIlate_dF3_BL DMIlate_dF3_C21 DMIlate_dF5_BL DMIlate_dF5_C21...
    Dprct1_BL Dprct1_C21 Dprct3_BL Dprct3_C21 Dprct5_BL Dprct5_C21 omean osum sem

%% Fig. 4d: Plot the % responsive PV and Pyr cells in Gq group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='Gq',:); 
means= [dataset.PercentResponsive(dataset.type=='PV',:)';...
    dataset.PercentResponsive(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PercentResponsive_SEM(dataset.type=='PV',:)';...
    dataset.PercentResponsive_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='PV',:);
dataset_C21_PV = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='PV',:);
dataset_BL_Pyr = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='Pyr',:);
dataset_C21_Pyr = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='Pyr',:);

plot([xdata(1,1),xdata(1,2)],[dataset_BL_PV.PercentResponsive,dataset_C21_PV.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot([xdata(2,1),xdata(2,2)],[dataset_BL_Pyr.PercentResponsive,dataset_C21_Pyr.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot(xdata(1,1),dataset_BL_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','[0.3 0.3 0.3]');
plot(xdata(1,2),dataset_C21_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
plot(xdata(2,1),dataset_BL_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','[0.3 0.3 0.3]');
plot(xdata(2,2),dataset_C21_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - %Responsive Cells'); % Choice of the name of the plot HERE
ylim([0 110]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.type=='PV',:);
normality = lillietest(dataset_PV.PercentResponsive) %0 normal
[h,p]=ttest(dataset_PV.PercentResponsive(dataset_PV.condition=='BL'),dataset_PV.PercentResponsive(dataset_PV.condition=='C21'))

dataset_Pyr=mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PercentResponsive) %0 normal
[h,p]=ttest(dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='BL'),dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='C21'))

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types
%% Fig. 4e-f: Plot the mean-evoked dF_F traces from all PV and Pyr cells in Gq group +/- C21
dataset = Data(Data.group=='Gq',:); 

BL_PV = []; C21_PV = []; BL_Pyr = []; C21_Pyr = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" && dataset.type(row)=='PV'
        BL_PV = [BL_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='PV'
        C21_PV = [C21_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "BL" && dataset.type(row)=='Pyr'
        BL_Pyr = [BL_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='Pyr'
        C21_Pyr = [C21_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL_PV = mean(cell2mat(BL_PV'), "omitnan"); std_BL_PV = std(cell2mat(BL_PV'), "omitnan")/sqrt(width(BL_PV)); % computation of the std error deviation
mean_C21_PV = mean(cell2mat(C21_PV'), "omitnan"); std_C21_PV = std(cell2mat(C21_PV'), "omitnan")/sqrt(width(C21_PV));
mean_BL_Pyr = mean(cell2mat(BL_Pyr'), "omitnan"); std_BL_Pyr = std(cell2mat(BL_Pyr'), "omitnan")/sqrt(width(BL_Pyr)); % computation of the std error deviation
mean_C21_Pyr = mean(cell2mat(C21_Pyr'), "omitnan"); std_C21_Pyr = std(cell2mat(C21_Pyr'), "omitnan")/sqrt(width(C21_Pyr));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

% Plot the traces with error bars for the standard errors - PV
means = [mean_BL_PV', mean_C21_PV']; stdErrors = [std_BL_PV', std_C21_PV'];
figure_meanEvokedResp_PV = figure;
axes_meanEvokedResp_PV = axes('Parent',figure_meanEvokedResp_PV);
hold(axes_meanEvokedResp_PV,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_PV,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_PV, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Mean whisker-evoked dF response of PV Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

% Plot the traces with error bars for the standard errors - Pyr
means = [mean_BL_Pyr', mean_C21_Pyr']; stdErrors = [std_BL_Pyr', std_C21_Pyr'];
figure_meanEvokedResp_Pyr = figure;
axes_meanEvokedResp_Pyr = axes('Parent',figure_meanEvokedResp_Pyr);
hold(axes_meanEvokedResp_Pyr,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_Pyr,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_Pyr, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Mean whisker-evoked dF response of Pyr Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

clear BL_PV BL_Pyr C21_PV C21_Pyr mean_BL_PV mean_BL_Pyr...
    mean_C21_PV mean_C21_Pyr std_BL_PV std_BL_Pyr std_C21_PV std_C21_Pyr...
    style means stdErrors x xLabels i dataset greybar

%% Fig. 4g: Plot the % responsive PV and Pyr cells in Gi group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='Gi',:); 
means= [dataset.PercentResponsive(dataset.type=='PV',:)';...
    dataset.PercentResponsive(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PercentResponsive_SEM(dataset.type=='PV',:)';...
    dataset.PercentResponsive_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='PV',:);
dataset_C21_PV = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='PV',:);
dataset_BL_Pyr = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='Pyr',:);
dataset_C21_Pyr = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='Pyr',:);

plot([xdata(1,1),xdata(1,2)],[dataset_BL_PV.PercentResponsive,dataset_C21_PV.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot([xdata(2,1),xdata(2,2)],[dataset_BL_Pyr.PercentResponsive,dataset_C21_Pyr.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot(xdata(1,1),dataset_BL_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
plot(xdata(1,2),dataset_C21_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
plot(xdata(2,1),dataset_BL_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
plot(xdata(2,2),dataset_C21_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - %Responsive Cells'); % Choice of the name of the plot HERE
ylim([0 110]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.type=='PV',:);
normality = lillietest(dataset_PV.PercentResponsive) %0 normal
[h,p]=ttest(dataset_PV.PercentResponsive(dataset_PV.condition=='BL'),dataset_PV.PercentResponsive(dataset_PV.condition=='C21'))

dataset_Pyr=mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PercentResponsive) %0 normal
[h,p]=ttest(dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='BL'),dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='C21'))

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types

%% Fig. 4h-i: Plot the mean-evoked dF_F traces from all PV and Pyr cells in Gi group +/- C21
dataset = Data(Data.group=='Gi',:); 

BL_PV = []; C21_PV = []; BL_Pyr = []; C21_Pyr = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" && dataset.type(row)=='PV'
        BL_PV = [BL_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='PV'
        C21_PV = [C21_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "BL" && dataset.type(row)=='Pyr'
        BL_Pyr = [BL_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='Pyr'
        C21_Pyr = [C21_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL_PV = mean(cell2mat(BL_PV'), "omitnan"); std_BL_PV = std(cell2mat(BL_PV'), "omitnan")/sqrt(width(BL_PV)); % computation of the std error deviation
mean_C21_PV = mean(cell2mat(C21_PV'), "omitnan"); std_C21_PV = std(cell2mat(C21_PV'), "omitnan")/sqrt(width(C21_PV));
mean_BL_Pyr = mean(cell2mat(BL_Pyr'), "omitnan"); std_BL_Pyr = std(cell2mat(BL_Pyr'), "omitnan")/sqrt(width(BL_Pyr)); % computation of the std error deviation
mean_C21_Pyr = mean(cell2mat(C21_Pyr'), "omitnan"); std_C21_Pyr = std(cell2mat(C21_Pyr'), "omitnan")/sqrt(width(C21_Pyr));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

% Plot the traces with error bars for the standard errors - PV
means = [mean_BL_PV', mean_C21_PV']; stdErrors = [std_BL_PV', std_C21_PV'];
figure_meanEvokedResp_PV = figure;
axes_meanEvokedResp_PV = axes('Parent',figure_meanEvokedResp_PV);
hold(axes_meanEvokedResp_PV,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_PV,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_PV, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Mean whisker-evoked dF response of PV Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

% Plot the traces with error bars for the standard errors - Pyr
means = [mean_BL_Pyr', mean_C21_Pyr']; stdErrors = [std_BL_Pyr', std_C21_Pyr'];
figure_meanEvokedResp_Pyr = figure;
axes_meanEvokedResp_Pyr = axes('Parent',figure_meanEvokedResp_Pyr);
hold(axes_meanEvokedResp_Pyr,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_Pyr,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_Pyr, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Mean whisker-evoked dF response of Pyr Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

clear BL_PV BL_Pyr C21_PV C21_Pyr mean_BL_PV mean_BL_Pyr...
    mean_C21_PV mean_C21_Pyr std_BL_PV std_BL_Pyr std_C21_PV std_C21_Pyr...
    style means stdErrors x xLabels i dataset greybar
%% Fig. 4j: Plot the % responsive PV and Pyr cells in mCherry group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='mCherry',:); 
means= [dataset.PercentResponsive(dataset.type=='PV',:)';...
    dataset.PercentResponsive(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PercentResponsive_SEM(dataset.type=='PV',:)';...
    dataset.PercentResponsive_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='PV',:);
dataset_C21_PV = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='PV',:);
dataset_BL_Pyr = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='BL' & mean_ByMouse.type=='Pyr',:);
dataset_C21_Pyr = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='C21' & mean_ByMouse.type=='Pyr',:);

plot([xdata(1,1),xdata(1,2)],[dataset_BL_PV.PercentResponsive,dataset_C21_PV.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot([xdata(2,1),xdata(2,2)],[dataset_BL_Pyr.PercentResponsive,dataset_C21_Pyr.PercentResponsive]',...
    '-ok','MarkerSize',12,'LineWidth',1,'Color','k','Marker','none');
plot(xdata(1,1),dataset_BL_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
plot(xdata(1,2),dataset_C21_PV.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
plot(xdata(2,1),dataset_BL_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
plot(xdata(2,2),dataset_C21_Pyr.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - %Responsive Cells'); % Choice of the name of the plot HERE
ylim([0 110]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.type=='PV',:);
normality = lillietest(dataset_PV.PercentResponsive) %0 normal
[h,p]=ttest(dataset_PV.PercentResponsive(dataset_PV.condition=='BL'),dataset_PV.PercentResponsive(dataset_PV.condition=='C21'))

dataset_Pyr=mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PercentResponsive) %0 normal
[h,p]=ttest(dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='BL'),dataset_Pyr.PercentResponsive(dataset_Pyr.condition=='C21'))

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types
%% Fig. S4k-l: Plot the mean-evoked dF_F traces from all PV and Pyr cells in mCherry group +/- C21
dataset = Data(Data.group=='mCherry',:); 

BL_PV = []; C21_PV = []; BL_Pyr = []; C21_Pyr = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" && dataset.type(row)=='PV'
        BL_PV = [BL_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='PV'
        C21_PV = [C21_PV,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "BL" && dataset.type(row)=='Pyr'
        BL_Pyr = [BL_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21"  && dataset.type(row)=='Pyr'
        C21_Pyr = [C21_Pyr,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL_PV = mean(cell2mat(BL_PV'), "omitnan"); std_BL_PV = std(cell2mat(BL_PV'), "omitnan")/sqrt(width(BL_PV)); % computation of the std error deviation
mean_C21_PV = mean(cell2mat(C21_PV'), "omitnan"); std_C21_PV = std(cell2mat(C21_PV'), "omitnan")/sqrt(width(C21_PV));
mean_BL_Pyr = mean(cell2mat(BL_Pyr'), "omitnan"); std_BL_Pyr = std(cell2mat(BL_Pyr'), "omitnan")/sqrt(width(BL_Pyr)); % computation of the std error deviation
mean_C21_Pyr = mean(cell2mat(C21_Pyr'), "omitnan"); std_C21_Pyr = std(cell2mat(C21_Pyr'), "omitnan")/sqrt(width(C21_Pyr));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)

% Plot the traces with error bars for the standard errors - PV
means = [mean_BL_PV', mean_C21_PV']; stdErrors = [std_BL_PV', std_C21_PV'];
figure_meanEvokedResp_PV = figure;
axes_meanEvokedResp_PV = axes('Parent',figure_meanEvokedResp_PV);
hold(axes_meanEvokedResp_PV,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_PV,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_PV, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Mean whisker-evoked dF response of PV Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

% Plot the traces with error bars for the standard errors - Pyr
means = [mean_BL_Pyr', mean_C21_Pyr']; stdErrors = [std_BL_Pyr', std_C21_Pyr'];
figure_meanEvokedResp_Pyr = figure;
axes_meanEvokedResp_Pyr = axes('Parent',figure_meanEvokedResp_Pyr);
hold(axes_meanEvokedResp_Pyr,'on');

style={'-k','-r'}; %Define colors for condition
greybar=NaN(size(mean_BL_Pyr,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp_Pyr, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Mean whisker-evoked dF response of Pyr Cells'); % Choice of the name of the plot HERE
ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])
legend('','','BL','','C21')

clear BL_PV BL_Pyr C21_PV C21_Pyr mean_BL_PV mean_BL_Pyr...
    mean_C21_PV mean_C21_Pyr std_BL_PV std_BL_Pyr std_C21_PV std_C21_Pyr...
    style means stdErrors x xLabels i dataset greybar
%% Fig. S9a: Plot baseline noise of all PV and Pyr cells in Gq group +/- C21

BL_PV = Data(Data.group=='Gq' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='Gq' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='Gq' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='Gq' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.Fbase_noise),mean(C21_PV.Fbase_noise);mean(BL_Pyr.Fbase_noise),mean(C21_Pyr.Fbase_noise)];
stdErrors = [std(BL_PV.Fbase_noise)/sqrt(size(BL_PV,1)),std(C21_PV.Fbase_noise)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.Fbase_noise)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.Fbase_noise)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Fbase_noise'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Noise'); % Choice of the name of the plot HERE
ylim([0 0.2]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='Gq' & Data.type=='PV',:);
normality = lillietest(dataset_PV.Fbase_noise) %1, not normal
lme = fitlme(dataset_PV,'Fbase_noise ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='Gq' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.Fbase_noise) %1, not normal
lme = fitlme(dataset_Pyr,'Fbase_noise ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h
%% Fig. S9b: Plot MI of all PV and Pyr cells in Gq group +/- C21

BL_PV = Data(Data.group=='Gq' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='Gq' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='Gq' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='Gq' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.MI_dF),mean(C21_PV.MI_dF);mean(BL_Pyr.MI_dF),mean(C21_Pyr.MI_dF)];
stdErrors = [std(BL_PV.MI_dF)/sqrt(size(BL_PV,1)),std(C21_PV.MI_dF)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.MI_dF)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.MI_dF)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on


scatter(xdata(1,1)*(ones(height(BL_PV),1)),BL_PV.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(height(C21_PV),1)),C21_PV.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(height(BL_Pyr),1)),BL_Pyr.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(height(C21_Pyr),1)),C21_Pyr.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('MI_dF'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - MI'); % Choice of the name of the plot HERE
ylim([-4 4.5]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='Gq' & Data.type=='PV',:);
normality = lillietest(dataset_PV.MI_dF) %1, not normal
lme = fitlme(dataset_PV,'MI_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='Gq' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.MI_dF) %1, not normal
lme = fitlme(dataset_Pyr,'MI_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h
%% Fig. S9c: Plot baseline noise of all PV and Pyr cells in Gi group +/- C21

BL_PV = Data(Data.group=='Gi' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='Gi' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='Gi' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='Gi' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.Fbase_noise),mean(C21_PV.Fbase_noise);mean(BL_Pyr.Fbase_noise),mean(C21_Pyr.Fbase_noise)];
stdErrors = [std(BL_PV.Fbase_noise)/sqrt(size(BL_PV,1)),std(C21_PV.Fbase_noise)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.Fbase_noise)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.Fbase_noise)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Fbase_noise'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Noise'); % Choice of the name of the plot HERE
ylim([0 0.2]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='Gi' & Data.type=='PV',:);
normality = lillietest(dataset_PV.Fbase_noise) %1, not normal
lme = fitlme(dataset_PV,'Fbase_noise ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='Gi' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.Fbase_noise) %1, not normal
lme = fitlme(dataset_Pyr,'Fbase_noise ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h
%% Fig. S9d: Plot MI of all PV and Pyr cells in Gi group +/- C21

BL_PV = Data(Data.group=='Gi' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='Gi' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='Gi' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='Gi' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.MI_dF),mean(C21_PV.MI_dF);mean(BL_Pyr.MI_dF),mean(C21_Pyr.MI_dF)];
stdErrors = [std(BL_PV.MI_dF)/sqrt(size(BL_PV,1)),std(C21_PV.MI_dF)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.MI_dF)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.MI_dF)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on


scatter(xdata(1,1)*(ones(height(BL_PV),1)),BL_PV.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(height(C21_PV),1)),C21_PV.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(height(BL_Pyr),1)),BL_Pyr.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(height(C21_Pyr),1)),C21_Pyr.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('MI_dF'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - MI'); % Choice of the name of the plot HERE
ylim([-4 4.5]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='Gi' & Data.type=='PV',:);
normality = lillietest(dataset_PV.MI_dF) %1, not normal
lme = fitlme(dataset_PV,'MI_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='Gi' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.MI_dF) %1, not normal
lme = fitlme(dataset_Pyr,'MI_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h
%% Fig. S9e: Plot early MI (0-1s post-stim) and late MI (1-2s post-stim) of all PV cells in Gi group +/- C21

BL_PV = Data(Data.group=='Gi' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='Gi' & Data.condition=='C21' & Data.type=='PV',:);

means= [mean(BL_PV.MIearly_dF),mean(C21_PV.MIearly_dF);mean(BL_PV.MIlate_dF),mean(C21_PV.MIlate_dF)];
stdErrors = [std(BL_PV.MIearly_dF)/sqrt(size(BL_PV,1)),std(C21_PV.MIearly_dF)/sqrt(size(C21_PV,1));...
    std(BL_PV.MIlate_dF)/sqrt(size(BL_PV,1)),std(C21_PV.MIlate_dF)/sqrt(size(C21_PV,1))];
conditions=unique(Data.condition);
xLabels={'early','late'};
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

scatter(xdata(1,1)*(ones(height(BL_PV),1)),BL_PV.MIearly_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(height(C21_PV),1)),C21_PV.MIearly_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(height(BL_PV),1)),BL_PV.MIlate_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(height(C21_PV),1)),C21_PV.MIlate_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('MI_dF'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - MIearly vs MIlate'); % Choice of the name of the plot HERE
ylim([-0.2 1]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics on the Early MI for PV cells
dataset_PV=Data(Data.group=='Gi' & Data.type=='PV',:);
normality = lillietest(dataset_PV.MIearly_dF) %1, not normal
lme = fitlme(dataset_PV,'MIearly_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics on the Late MI for PV cells
dataset_PV=Data(Data.group=='Gi' & Data.type=='PV',:);
normality = lillietest(dataset_PV.MIlate_dF) %1, not normal
lme = fitlme(dataset_PV,'MIlate_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata h lme
%% Fig. S9f: Plot baseline noise of all PV and Pyr cells in mCherry group +/- C21

BL_PV = Data(Data.group=='mCherry' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='mCherry' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='mCherry' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='mCherry' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.Fbase_noise),mean(C21_PV.Fbase_noise);mean(BL_Pyr.Fbase_noise),mean(C21_Pyr.Fbase_noise)];
stdErrors = [std(BL_PV.Fbase_noise)/sqrt(size(BL_PV,1)),std(C21_PV.Fbase_noise)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.Fbase_noise)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.Fbase_noise)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Fbase_noise'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Noise'); % Choice of the name of the plot HERE
ylim([0 0.2]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='mCherry' & Data.type=='PV',:);
normality = lillietest(dataset_PV.Fbase_noise) %1, not normal
lme = fitlme(dataset_PV,'Fbase_noise ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='mCherry' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.Fbase_noise) %1, not normal
lme = fitlme(dataset_Pyr,'Fbase_noise ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h

%% Fig S9g: Plot MI of all PV and Pyr cells in mCherry group +/- C21

BL_PV = Data(Data.group=='mCherry' & Data.condition=='BL' & Data.type=='PV',:);
C21_PV = Data(Data.group=='mCherry' & Data.condition=='C21' & Data.type=='PV',:);
BL_Pyr = Data(Data.group=='mCherry' & Data.condition=='BL' & Data.type=='Pyr',:);
C21_Pyr = Data(Data.group=='mCherry' & Data.condition=='C21' & Data.type=='Pyr',:);

means= [mean(BL_PV.MI_dF),mean(C21_PV.MI_dF);mean(BL_Pyr.MI_dF),mean(C21_Pyr.MI_dF)];
stdErrors = [std(BL_PV.MI_dF)/sqrt(size(BL_PV,1)),std(C21_PV.MI_dF)/sqrt(size(C21_PV,1));...
    std(BL_Pyr.MI_dF)/sqrt(size(BL_Pyr,1)),std(C21_Pyr.MI_dF)/sqrt(size(C21_Pyr,1))];
conditions=unique(Data.condition);
types=unique(Data.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on


scatter(xdata(1,1)*(ones(height(BL_PV),1)),BL_PV.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(height(C21_PV),1)),C21_PV.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(height(BL_Pyr),1)),BL_Pyr.MI_dF,...
    'ok','MarkerFaceColor',[0.3 0.3 0.3],'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(height(C21_Pyr),1)),C21_Pyr.MI_dF,...
    'ok','MarkerFaceColor','r','XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('MI_dF'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - MI'); % Choice of the name of the plot HERE
ylim([-4 4.5]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')
pbaspect([1 1.5 1])

% Get statistics for PV cells
dataset_PV=Data(Data.group=='mCherry' & Data.type=='PV',:);
normality = lillietest(dataset_PV.MI_dF) %1, not normal
lme = fitlme(dataset_PV,'MI_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

% Get statistics for Pyr cells
dataset_Pyr=Data(Data.group=='mCherry' & Data.type=='Pyr',:);
normality = lillietest(dataset_Pyr.MI_dF) %1, not normal
lme = fitlme(dataset_Pyr,'MI_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear means stdErrors conditions xLabels colors BL_PV C21_PV dataset_PV...
    dataset_Pyr normality BL_Pyr C21_Pyr xdata types lme h

%% Fig. S10a,d,g: Plot change in % whisker responsive PV and Pyr Cells for Gi, Gq, and mCherry
means= [delta_ByGroup1.Dprct(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gi');...
    delta_ByGroup1.Dprct(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gi')];
stdErrors= [delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gi');...
    delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gi')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='Gi' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='Gi' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='Gi' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='Gi' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='Gi' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='Gi' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in %Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Delta %Responsive Cells'); % Choice of the name of the plot HERE
ylim([-50 50]) % add a limit here if needed to scale different plots together reliability
legend(doses)

%Plot for Gq
means= [delta_ByGroup1.Dprct(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gq');...
    delta_ByGroup1.Dprct(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gq')];
stdErrors= [delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gq');...
    delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gq')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

figure
h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='Gq' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='Gq' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='Gq' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='Gq' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='Gq' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='Gq' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in %Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Delta %Responsive Cells'); % Choice of the name of the plot HERE
ylim([-100 50]) % add a limit here if needed to scale different plots together reliability
legend(doses)

%Plot for mCherry
means= [delta_ByGroup1.Dprct(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='mCherry');...
    delta_ByGroup1.Dprct(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.Dprct(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.Dprct(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='mCherry')];
stdErrors= [delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='mCherry');...
    delta_ByGroup1.Dprct_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.Dprct_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.Dprct_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='mCherry')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

figure
h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='mCherry' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='mCherry' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='mCherry' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.Dprct(delta_ByMouse1.group=='mCherry' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.Dprct(delta_ByMouse3.group=='mCherry' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.Dprct(delta_ByMouse5.group=='mCherry' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in %Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Delta %Responsive Cells'); % Choice of the name of the plot HERE
ylim([-60 50]) % add a limit here if needed to scale different plots together reliability
legend(doses)

clear means stdErrors types xdata xLabels doses h i
%% Fig. S10b,e,h: Plot change in MI for PV and Pyr cells for Gq, Gi, and mCherry pre/post C21
%Plot for Gi
%Use MIearly_dF for PV and MI_dF for Pyr to match main figure
means= [delta_ByGroup1.DMIearly_dF(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.DMIearly_dF(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.DMIearly_dF(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gi');...
    delta_ByGroup1.DMI_dF(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.DMI_dF(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.DMI_dF(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gi')];
stdErrors= [delta_ByGroup1.DMIearly_dF_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.DMIearly_dF_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.DMIearly_dF_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gi');...
    delta_ByGroup1.DMI_dF_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gi'),...
    delta_ByGroup3.DMI_dF_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gi'),...
    delta_ByGroup5.DMI_dF_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gi')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.DMIearly_dF(delta_ByMouse1.group=='Gi' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.DMIearly_dF(delta_ByMouse3.group=='Gi' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.DMIearly_dF(delta_ByMouse5.group=='Gi' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.DMI_dF(delta_ByMouse1.group=='Gi' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.DMI_dF(delta_ByMouse3.group=='Gi' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.DMI_dF(delta_ByMouse5.group=='Gi' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in MI_dF Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Delta MI_dF all Cells'); % Choice of the name of the plot HERE
ylim([-0.25 0.35]) % add a limit here if needed to scale different plots together reliability
legend(doses)


%Plot for Gq
means= [delta_ByGroup1.DMI_dF(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.DMI_dF(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.DMI_dF(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gq');...
    delta_ByGroup1.DMI_dF(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.DMI_dF(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.DMI_dF(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gq')];
stdErrors= [delta_ByGroup1.DMI_dF_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.DMI_dF_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.DMI_dF_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='Gq');...
    delta_ByGroup1.DMI_dF_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='Gq'),...
    delta_ByGroup3.DMI_dF_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='Gq'),...
    delta_ByGroup5.DMI_dF_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='Gq')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

figure
h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.DMI_dF(delta_ByMouse1.group=='Gq' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.DMI_dF(delta_ByMouse3.group=='Gq' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.DMI_dF(delta_ByMouse5.group=='Gq' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.DMI_dF(delta_ByMouse1.group=='Gq' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.DMI_dF(delta_ByMouse3.group=='Gq' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.DMI_dF(delta_ByMouse5.group=='Gq' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in MI_dF Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Delta MI_dF all Cells'); % Choice of the name of the plot HERE
ylim([-0.25 0.35]) % add a limit here if needed to scale different plots together reliability
legend(doses)

%Plot for mCherry
%Use MIearly_dF for PV and MI_dF for Pyr to match main figure
means= [delta_ByGroup1.DMI_dF(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.DMI_dF(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.DMI_dF(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='mCherry');...
    delta_ByGroup1.DMI_dF(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.DMI_dF(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.DMI_dF(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='mCherry')];
stdErrors= [delta_ByGroup1.DMI_dF_SEM(delta_ByGroup1.type=='PV' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.DMI_dF_SEM(delta_ByGroup3.type=='PV' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.DMI_dF_SEM(delta_ByGroup5.type=='PV' & delta_ByGroup5.group=='mCherry');...
    delta_ByGroup1.DMI_dF_SEM(delta_ByGroup1.type=='Pyr' & delta_ByGroup1.group=='mCherry'),...
    delta_ByGroup3.DMI_dF_SEM(delta_ByGroup3.type=='Pyr' & delta_ByGroup3.group=='mCherry'),...
    delta_ByGroup5.DMI_dF_SEM(delta_ByGroup5.type=='Pyr' & delta_ByGroup5.group=='mCherry')];
doses={'1 mg/kg','3 mg/kg', '5 mg/kg'};
types=unique(delta_ByGroup1.type);
xLabels=types;

figure
h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[1 0.9 0.9];
h(2).FaceColor=[1 0.6 0.6];
h(3).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints;h(3).XEndPoints]';
hold on

plot(xdata(1,1),delta_ByMouse1.DMI_dF(delta_ByMouse1.group=='mCherry' & delta_ByMouse1.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(1,2),delta_ByMouse3.DMI_dF(delta_ByMouse3.group=='mCherry' & delta_ByMouse3.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(1,3),delta_ByMouse5.DMI_dF(delta_ByMouse5.group=='mCherry' & delta_ByMouse5.type=='PV',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);
plot(xdata(2,1),delta_ByMouse1.DMI_dF(delta_ByMouse1.group=='mCherry' & delta_ByMouse1.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.9 0.9]);
plot(xdata(2,2),delta_ByMouse3.DMI_dF(delta_ByMouse3.group=='mCherry' & delta_ByMouse3.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0.6 0.6]);
plot(xdata(2,3),delta_ByMouse5.DMI_dF(delta_ByMouse5.group=='mCherry' & delta_ByMouse5.type=='Pyr',:),...
    'ok','MarkerSize',12,'MarkerFaceColor',[1 0 0]);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2,'CapSize',10);
errorbar(xdata(:,3),means(:,3),stdErrors(:,3),'k','linestyle','none','LineWidth',2,'CapSize',10);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Change in MI_dF Cells'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Delta MI_dF all Cells'); % Choice of the name of the plot HERE
ylim([-0.25 0.35]) % add a limit here if needed to scale different plots together reliability
legend(doses)

clear means stdErrors types xdata xLabels doses

%% Section 6 (Required): Load the spontaneous data for PV DREADD Experiments - Analysis of PV and Pyr cells
clear %clear workspace

% Load the data
load('/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/Acute_DREADDs_Suite2P/C21_5mg_Kg/Spontaneous Data Analysis/Acute_DREADDs_5_SpontDataByCell_spks.mat')

% Calculate peaks for dF_F
framerate=10.0143333333333;
Data=Spont_Data_ByCell;

for ii=1:height(Data)
    %Calculate the baseline of dF_F
    Fraw=double(Data.F{ii});
    Fneu=double(Data.Fneu{ii});
    F=(Fraw-0.7*Fneu)';
    Fsize=size(F);
    baselineframes=round(10*framerate);

    %-----------Calculate Baseline F------------------------
    basecalc=zeros(1,Fsize(1)-baselineframes+1);
    for i=1:Fsize(1)-baselineframes+1 %find 10 seconds of baseline. RECALC based on framerate
        basecalc(i)=std(F(i:i+baselineframes-1)); %checks STD for all values
    end
    [~,mindex]=min(basecalc); %calculate the position of lowest deviation for baseline
    dF_std=std(Data.dF_F{ii}(mindex:mindex+baselineframes-1));
    Data.dF_std(ii)=dF_std;
    %----------------End Calculate Baseline F--------------------------------
    clear Fraw Fneu F Fsize baselineframes basecalc mindex
    
    %Calculate the mean Z-score for the entire trace
    Data.mean_dF(ii)=mean(Data.dF_F{ii});
    
    %Calculate the AUC of the entire trace
    Data.AUC_dF(ii)=trapz(Data.dF_F{ii})/framerate;

    %Peak finding for entire trace
    warning('off','signal:findpeaks:largeMinPeakHeight')
    
    %Set MinPeakHeight to 4x STD of baseline and MinPeakProminence to 3x
    %STD to match the Z-scored analysis settings
    [pks,locs,w,p]=findpeaks(Data.dF_F{ii},'MinPeakHeight',(4*dF_std),'MinPeakDistance',round(framerate*0.5),'MinPeakProminence',(3*dF_std));
%     [pks,locs,w,p]=findpeaks(Data.dF_F{ii},'MinPeakDistance',round(framerate*0.5));
    Data.Peaks_Amp_dF{ii}=pks;
    Data.Peaks_Loc_dF{ii}=locs;
    Data.Peaks_Width_dF{ii}=w;
    Data.Peaks_Prom_dF{ii}=p;

    traceduration=size(Data.dF_F{ii},2)/framerate;
    
    Data.PeaksAmp_mean_dF(ii)=mean(pks); %mean peak amplitude Z-score
    Data.PeakFreq_dF(ii)=size(locs,2)/traceduration; %peak frequency as peaks/sec
    clear pks locs w p traceduration dF_std
end

warning('on','signal:findpeaks:largeMinPeakHeight')
Spont_Data_ByCell=Data;

% Pool all cells by group
T1=Spont_Data_ByCell(:,{'condition','group','type'});
[G,mean_allcells]=findgroups(T1);
C=groupcounts(T1,{'condition','group','type'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_allcells.TotalCells=C.GroupCount;
mean_allcells.AUC=splitapply(omean,Spont_Data_ByCell.AUC,G);
mean_allcells.AUC_SEM=splitapply(sem,Spont_Data_ByCell.AUC,G);
mean_allcells.mean_Z=splitapply(omean,Spont_Data_ByCell.mean_Z,G);
mean_allcells.mean_Z_SEM=splitapply(sem,Spont_Data_ByCell.mean_Z,G);
mean_allcells.PeaksAmp_mean=splitapply(omean,Spont_Data_ByCell.PeaksAmp_mean,G);
mean_allcells.PeaksAmp_mean_SEM=splitapply(sem,Spont_Data_ByCell.PeaksAmp_mean,G);
mean_allcells.PeakFreq=splitapply(omean,Spont_Data_ByCell.PeakFreq,G);
mean_allcells.PeakFreq_SEM=splitapply(sem,Spont_Data_ByCell.PeakFreq,G);

% For dF_F
mean_allcells.AUC_dF=splitapply(omean,Spont_Data_ByCell.AUC_dF,G);
mean_allcells.AUC_dF_SEM=splitapply(sem,Spont_Data_ByCell.AUC_dF,G);
mean_allcells.mean_dF=splitapply(omean,Spont_Data_ByCell.mean_dF,G);
mean_allcells.mean_dF_SEM=splitapply(sem,Spont_Data_ByCell.mean_dF,G);
mean_allcells.PeaksAmp_mean_dF=splitapply(omean,Spont_Data_ByCell.PeaksAmp_mean_dF,G);
mean_allcells.PeaksAmp_mean_dF_SEM=splitapply(sem,Spont_Data_ByCell.PeaksAmp_mean_dF,G);
mean_allcells.PeakFreq_dF=splitapply(omean,Spont_Data_ByCell.PeakFreq_dF,G);
mean_allcells.PeakFreq_dF_SEM=splitapply(sem,Spont_Data_ByCell.PeakFreq_dF,G);

clear C Data framerate G i ii omean omean2 sem sem2 T1
%% Fig. S11a: Plot the AUC of the spontaneous activity trace in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= [dataset.AUC_dF(dataset.type=='PV',:)';...
    dataset.AUC_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.AUC_dF_SEM(dataset.type=='PV',:)';...
    dataset.AUC_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('AUC (dF*sec)'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - AUC_dF'); % Choice of the name of the plot HERE
% ylim([-600 1500]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.AUC_dF) %0 normal
lme = fitlme(dataset_PV,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.AUC_dF) %1 not normal
lme = fitlme(dataset_Pyr,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme
%% Fig. S11b: Plot the Peak Frequency of spontaneous activity in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= [dataset.PeakFreq_dF(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeakFreq_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Frequency (Hz)'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - PeakFreq_dF'); % Choice of the name of the plot HERE
% ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeakFreq_dF) %0 normal
lme = fitlme(dataset_PV,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme

%% Fig. S11c: Plot the mean Peak Amplitude of spontaneous activity in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= [dataset.PeaksAmp_mean_dF(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeaksAmp_mean_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Amplitude (dF_F)'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - PeaksAmp_mean_dF'); % Choice of the name of the plot HERE
% ylim([0 30]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeaksAmp_mean_dF) %0 normal
lme = fitlme(dataset_PV,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme
%% Fig. S11d: Plot the AUC of the spontaneous activity trace in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:); 
means= [dataset.AUC_dF(dataset.type=='PV',:)';...
    dataset.AUC_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.AUC_dF_SEM(dataset.type=='PV',:)';...
    dataset.AUC_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('AUC (dF_F*sec)'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - AUC_dF'); % Choice of the name of the plot HERE
% ylim([-400 1400]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.AUC_dF) %0 normal
lme = fitlme(dataset_PV,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.AUC_dF) %1 not normal
lme = fitlme(dataset_Pyr,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme
%% Fig. S11e: Plot the Peak Frequency of spontaneous activity in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:); 
means= [dataset.PeakFreq_dF(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeakFreq_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Frequency (Hz)'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - PeakFreq_dF'); % Choice of the name of the plot HERE
% ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeakFreq_dF) %1, not normal
lme = fitlme(dataset_PV,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme

%% Fig. S11f: Plot the mean Peak Amplitude of spontaneous activity in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:); 
means= [dataset.PeaksAmp_mean_dF(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeaksAmp_mean_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Amplitude (dF_F)'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - PeaksAmp_mean'); % Choice of the name of the plot HERE
% ylim([0 30]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeaksAmp_mean_dF) %1, not normal
lme = fitlme(dataset_PV,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme
%% Fig. S11g: Plot the AUC of the spontaneous activity trace in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= [dataset.AUC_dF(dataset.type=='PV',:)';...
    dataset.AUC_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.AUC_dF_SEM(dataset.type=='PV',:)';...
    dataset.AUC_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('AUC_dF (Z-score*sec)'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - AUC_dF'); % Choice of the name of the plot HERE
% ylim([-400 1200]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.AUC_dF) %0 normal
lme = fitlme(dataset_PV,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.AUC_dF) %1 not normal
lme = fitlme(dataset_Pyr,'AUC_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme

%% Fig. S11h: Plot the Peak Frequency of spontaneous activity in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= [dataset.PeakFreq_dF(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeakFreq_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeakFreq_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Frequency (Hz)'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - PeakFreq_dF'); % Choice of the name of the plot HERE
% ylim([0 0.8]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeakFreq_dF) %0 normal
lme = fitlme(dataset_PV,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeakFreq_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme

%% Fig. S11i: Plot the mean Peak Amplitude of spontaneous activity in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= [dataset.PeaksAmp_mean_dF(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF(dataset.type=='Pyr',:)'];
stdErrors = [dataset.PeaksAmp_mean_dF_SEM(dataset.type=='PV',:)';...
    dataset.PeaksAmp_mean_dF_SEM(dataset.type=='Pyr',:)'];
conditions=unique(dataset.condition);
types=unique(dataset.type);
xLabels=types;
colors=[0.3 0.3 0.3; 1 0 0];

h=bar(means,'FaceColor','flat','FaceAlpha',0.5)
h(1).FaceColor=[0.3 0.3 0.3];
h(2).FaceColor=[1 0 0];
xdata=[h(1).XEndPoints;h(2).XEndPoints]';
hold on

dataset_BL_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='PV',:);
dataset_C21_PV = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='PV',:);
dataset_BL_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL' & Spont_Data_ByCell.type=='Pyr',:);
dataset_C21_Pyr = Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21' & Spont_Data_ByCell.type=='Pyr',:);

scatter(xdata(1,1)*(ones(size(dataset_BL_PV,1),1)),...
    dataset_BL_PV.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(1,2)*(ones(size(dataset_C21_PV,1),1)),...
    dataset_C21_PV.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,1)*(ones(size(dataset_BL_Pyr,1),1)),...
    dataset_BL_Pyr.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(xdata(2,2)*(ones(size(dataset_C21_Pyr,1),1)),...
    dataset_C21_Pyr.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(xdata(:,1),means(:,1),stdErrors(:,1),'k','linestyle','none','LineWidth',2);
errorbar(xdata(:,2),means(:,2),stdErrors(:,2),'k','linestyle','none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Cell Type'); 
ylabel('Peak Amplitude (dF_F)'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - PeaksAmp_mean_dF'); % Choice of the name of the plot HERE
% ylim([0 30]) % add a limit here if needed to scale different plots together reliability
legend('BL','C21')

%Get statistics
dataset_PV=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='PV',:);
normality = lillietest(dataset_PV.PeaksAmp_mean_dF) %0 normal
lme = fitlme(dataset_PV,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

dataset_Pyr=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.type=='Pyr',:);
normality = lillietest(dataset_Pyr.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_Pyr,'PeaksAmp_mean_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL_PV dataset_C21_PV dataset_PV...
    dataset_Pyr normality h p dataset_BL_Pyr dataset_C21_Pyr xdata types lme
%% Section 7 (Required): Load the whisker evoked data for PV DREADD Experiments - Analysis of SST cells
clear %clear workspace variables

load("/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/SST_DREADDs_Cre-Flip/Spike Deconvolution Analysis/SST_DREADDs_5_DataByCell_SpikeProb.mat")

% Get mean evoked dF_F traces
Data=Data_ByCell;
allstimsallROIs_dF=cell(size(Data,1),1);
framerate=10.0143333333333;
for k=1:size(Data,1)
    stimulus_trace=Data.stim_epochs{k};
    signal_traces=Data.dF_F{k};
    
    idx=find(stimulus_trace);
    stimstarts=idx(1,1);
    stimends=[];
    for i=2:size(idx,1)
        if idx(i)>(idx(i-1)+1)
            stimstarts=[stimstarts,idx(i)];
            stimends=[stimends,idx(i-1)];
        else
        end
    end
    stimends=[stimends,idx(end)];
    stimstarts=stimstarts';
    stimends=stimends';

    numberofstims=size(stimstarts,1);
    timeafterstim = 3.1; %Given an inter-stim interval of 3.5-6 seconds 
    framesafterstim = round(timeafterstim*framerate);
    timebeforestim = 1; %To be used as a pre-stimulus baseline
    framesbeforestim = round(timebeforestim*framerate);
    chopsize = framesbeforestim+framesafterstim+1;
    peakcut=framesbeforestim+round(2*framerate); %Set to baseline plus 2 seconds

    ROIcount = size(signal_traces,1); 

    %creates a 3D matrix of all stimulus-evoked epochs. 3 parameters: ROIs, stimulation (from
    %1 to numberofstims), frames (deltaF data)
    allstimsallROIs = zeros(chopsize,numberofstims);

    %Now indexing into this 3D matrix. Frames x Stims x ROI
    for j = 1:numberofstims %establishes the matrix columns, i.e. the whisker stimulations
        allstimsallROIs(:,j) = (signal_traces((stimstarts(j)-framesbeforestim):(stimstarts(j)+framesafterstim)));
    end  
    allstimsallROIs_dF{k}=allstimsallROIs;
end

for i=1:size(allstimsallROIs_dF,1)
    allstimsallROIs_dF_mean{i}=mean(allstimsallROIs_dF{i},2)';
end
allstimsallROIs_dF_mean=allstimsallROIs_dF_mean';

%Add to the table
Data.allstimsallROIs_dF=allstimsallROIs_dF;
Data.allstimsallROIs_dF_mean=allstimsallROIs_dF_mean;

%Calculate MI using dF
for i=1:height(Data)
    mean_trace=Data.allstimsallROIs_dF_mean{i};
    Data.MI_dF(i)=mean(mean_trace(11:30))-mean(mean_trace(1:10));
    Data.MIearly_dF(i)=mean(mean_trace(11:20))-mean(mean_trace(1:10));
    Data.MIlate_dF(i)=mean(mean_trace(21:30))-mean(mean_trace(1:10));
end

Data_ByCell=Data;

clear allstimsallROIs_dF allstimsallROIs_dF_mean stimulus_trace signal_traces...
    stimends stimstarts timeafterstim timebeforestim peakcut i j k framesbeforestim...
    framesafterstim chopsize idx ROIcount

% Get means for each condition/group from individual mice
Data=Data_ByCell;
omean=@(x) mean(x,'omitnan');
osum=@(x) sum(x);
T1=Data(:,{'MouseID','condition','group'});
[G,mean_ByMouse]=findgroups(T1);
C=groupcounts(T1,{'MouseID','condition','group'});

mean_ByMouse.NumRespCells=splitapply(osum,Data.Responsivity,G);
mean_ByMouse.TotalCells=C.GroupCount;
mean_ByMouse.PercentResponsive=(mean_ByMouse.NumRespCells./mean_ByMouse.TotalCells)*100;

%Then summarize by group to get mean with n=# of mice
omean=@(x) mean(x,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
T1=mean_ByMouse(:,{'condition','group'});
[G,mean_ByGroup]=findgroups(T1);

mean_ByGroup.PercentResponsive=splitapply(omean,mean_ByMouse.PercentResponsive,G);
mean_ByGroup.PercentResponsive_SEM=splitapply(sem,mean_ByMouse.PercentResponsive,G);

% Pool all cells by group for responsive cells and non-responsive cells
Data_Resp=Data(Data.Responsivity==1,:);
T1=Data_Resp(:,{'condition','group'});
[G,mean_Resp]=findgroups(T1);
C=groupcounts(T1,{'condition','group'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_Resp.TotalCells=C.GroupCount;
mean_Resp.AUC=splitapply(omean,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.AUC_SEM=splitapply(sem,Data_Resp.allstimsallROIs_AUC,G);
mean_Resp.MI=splitapply(omean,Data_Resp.MI,G);
mean_Resp.MI_SEM=splitapply(sem,Data_Resp.MI,G);
mean_Resp.Amp=splitapply(omean,Data_Resp.PeaksMax_mean,G);
mean_Resp.Amp_SEM=splitapply(sem,Data_Resp.PeaksMax_mean,G);
mean_Resp.Latency=splitapply(omean,Data_Resp.PeaksLatency_mean,G);
mean_Resp.Latency_SEM=splitapply(sem,Data_Resp.PeaksLatency_mean,G);
mean_Resp.PrctStimActive=splitapply(omean,Data_Resp.PrctStimActive,G);
mean_Resp.PrctStimActive_SEM=splitapply(sem,Data_Resp.PrctStimActive,G);
mean_Resp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_Resp.allstimsROIs_mean),G);
mean_Resp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_Resp.allstimsROIs_mean),G);

Data_NResp=Data(Data.Responsivity==0,:);
T1=Data_NResp(:,{'condition','group'});
[G,mean_NResp]=findgroups(T1);
C=groupcounts(T1,{'condition','group'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_NResp.TotalCells=C.GroupCount;
mean_NResp.AUC=splitapply(omean,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.AUC_SEM=splitapply(sem,Data_NResp.allstimsallROIs_AUC,G);
mean_NResp.MI=splitapply(omean,Data_NResp.MI,G);
mean_NResp.MI_SEM=splitapply(sem,Data_NResp.MI,G);
mean_NResp.Amp=splitapply(omean,Data_NResp.PeaksMax_mean,G);
mean_NResp.Amp_SEM=splitapply(sem,Data_NResp.PeaksMax_mean,G);
mean_NResp.Latency=splitapply(omean,Data_NResp.PeaksLatency_mean,G);
mean_NResp.Latency_SEM=splitapply(sem,Data_NResp.PeaksLatency_mean,G);
mean_NResp.PrctStimActive=splitapply(omean,Data_NResp.PrctStimActive,G);
mean_NResp.PrctStimActive_SEM=splitapply(sem,Data_NResp.PrctStimActive,G);
mean_NResp.allstimsROIs_mean=splitapply(omean2,cell2mat(Data_NResp.allstimsROIs_mean),G);
mean_NResp.allstimsROIs_SEM=splitapply(sem2,cell2mat(Data_NResp.allstimsROIs_mean),G);

clear allstimsallROIs C framerate G mean_trace numberofstims omean...
    omean2 osum sem sem2 T1
%% Fig. 5c: Plot the % responsive SST cells in Gq group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='Gq',:); 
means= dataset.PercentResponsive;
stdErrors = dataset.PercentResponsive_SEM;
conditions=unique(dataset.condition);
xLabels=conditions;
colors=[0.3 0.3 0.3; 1 0 0];

figure_PrctResp_Gq=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp_Gq.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',1); % parameters to modify for the design (if needed)

dataset_BL = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='BL',:);
dataset_C21 = mean_ByMouse(mean_ByMouse.group=='Gq' & mean_ByMouse.condition=='C21',:);

plot([dataset_BL.PercentResponsive,dataset_C21.PercentResponsive]','-ok','MarkerSize',12,...
    'LineWidth',1,'Color','k','Marker','none');
plot(1,dataset_BL.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','k');
plot(2,dataset_C21.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Condition'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - %Responsive SSTs'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset2=mean_ByMouse(mean_ByMouse.group=='Gq',:);
normality = lillietest(dataset2.PercentResponsive) %1 not normal
[p,h]=signrank(dataset2.PercentResponsive(dataset2.condition=='BL'),dataset2.PercentResponsive(dataset2.condition=='C21'))

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21 dataset2...
    normality h p k

%% Fig. 5d: Plot the mean-evoked dF_F traces from all SST cells in Gq group
dataset = Data_ByCell(Data_ByCell.group=='Gq',:); 

BL = []; C21 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" 
        BL = [BL,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21" 
        C21 = [C21,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL = mean(cell2mat(BL'), "omitnan"); std_BL = std(cell2mat(BL'), "omitnan")/sqrt(width(BL)); % computation of the std error deviation
mean_C21 = mean(cell2mat(C21'), "omitnan"); std_C21 = std(cell2mat(C21'), "omitnan")/sqrt(width(C21));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_BL', mean_C21']; stdErrors = [std_BL', std_C21'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-k','-r'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_BL,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM3DGq - Mean whisker-evoked dF_F response of all cells SSTs'); % Choice of the name of the plot HERE
ylim([0 0.14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

%Get statistics
dataset = Data_ByCell(Data_ByCell.group=='Gq',:); 
dataset.condition = reordercats(dataset.condition, {'BL', 'C21'}); 
normality = lillietest(dataset.MI) %1, not normal
lme = fitlme(dataset,'MI ~ condition + (1 | MouseID)','FitMethod','REML')


clear BL C21 mean_BL mean_C21 std_BL std_C21 means stdErrors x xLabels greybar...
    i style dataset normality lme
%% Fig. 5e: Plot the % responsive SST cells in Gi group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='Gi',:); 
means= dataset.PercentResponsive;
stdErrors = dataset.PercentResponsive_SEM;
conditions=unique(dataset.condition);
xLabels=conditions;
colors=[0.3 0.3 0.3; 1 0 0];

figure_PrctResp=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',1); % parameters to modify for the design (if needed)

dataset_BL = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='BL',:);
dataset_C21 = mean_ByMouse(mean_ByMouse.group=='Gi' & mean_ByMouse.condition=='C21',:);

plot([dataset_BL.PercentResponsive,dataset_C21.PercentResponsive]','-ok','MarkerSize',12,...
    'LineWidth',1,'Color','k','Marker','none');
plot(1,dataset_BL.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','k');
plot(2,dataset_C21.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Condition'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - %Responsive SSTs'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset2=mean_ByMouse(mean_ByMouse.group=='Gi',:);
normality = lillietest(dataset2.PercentResponsive) %1 not normal
[p,h]=signrank(dataset2.PercentResponsive(dataset2.condition=='BL'),dataset2.PercentResponsive(dataset2.condition=='C21'));

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21 dataset2...
    normality h p k

%% Fig. 5f: Plot the mean-evoked dF_F traces from all SST cells in Gi group +/- C21
dataset = Data_ByCell(Data_ByCell.group=='Gi',:); 

BL = []; C21 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" 
        BL = [BL,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21" 
        C21 = [C21,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL = mean(cell2mat(BL'), "omitnan"); std_BL = std(cell2mat(BL'), "omitnan")/sqrt(width(BL)); % computation of the std error deviation
mean_C21 = mean(cell2mat(C21'), "omitnan"); std_C21 = std(cell2mat(C21'), "omitnan")/sqrt(width(C21));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_BL', mean_C21']; stdErrors = [std_BL', std_C21'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-k','-r'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_BL,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-HM4DGi - Mean whisker-evoked dF_F response of all cells SSTs'); % Choice of the name of the plot HERE
ylim([0 0.14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

%Get statistics
dataset = Data_ByCell(Data_ByCell.group=='Gi',:); 
dataset.condition = reordercats(dataset.condition, {'BL', 'C21'}); 
normality = lillietest(dataset.MIearly_dF) %1, not normal
lme = fitlme(dataset,'MIearly_dF ~ condition + (1 | MouseID)','FitMethod','REML')

normality = lillietest(dataset.MIlate_dF) %1, not normal
lme = fitlme(dataset,'MIlate_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear BL C21 mean_BL mean_C21 std_BL std_C21 means stdErrors x xLabels dataset...
    greybar i style normality lme
%% Fig. 5g: Plot the % responsive SST cells in mCherry group +/- C21
dataset = mean_ByGroup(mean_ByGroup.group=='mCherry',:); 
means= dataset.PercentResponsive;
stdErrors = dataset.PercentResponsive_SEM;
conditions=unique(dataset.condition);
xLabels=conditions;
colors=[0.3 0.3 0.3; 1 0 0];

figure_PrctResp_mC=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PrctResp_mC.CData(k,:) = colors(k,:);
end
errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',1); % parameters to modify for the design (if needed)

dataset_BL = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='BL',:);
dataset_C21 = mean_ByMouse(mean_ByMouse.group=='mCherry' & mean_ByMouse.condition=='C21',:);

plot([dataset_BL.PercentResponsive,dataset_C21.PercentResponsive]','-ok','MarkerSize',12,...
    'LineWidth',1,'Color','k','Marker','none');
plot(1,dataset_BL.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','k');
plot(2,dataset_C21.PercentResponsive,...
    'ok','MarkerSize',12,'MarkerFaceColor','r');
hold off;

% Add the x-axis labels and title
set(gca, 'XTick', 1:length(xLabels), 'XTickLabel', xLabels);
xlabel('Condition'); 
ylabel('% Responsive Cells'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - %Responsive SSTs'); % Choice of the name of the plot HERE
ylim([0 100]) % add a limit here if needed to scale different plots together reliability
pbaspect([1 2.13 1])

%Get statistics
dataset2=mean_ByMouse(mean_ByMouse.group=='mCherry',:);
normality = lillietest(dataset2.PercentResponsive) %0 normal
[h,p]=ttest(dataset2.PercentResponsive(dataset2.condition=='BL'),dataset2.PercentResponsive(dataset2.condition=='C21'))

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21 dataset2...
    normality h p k

%% Fig. 5h: Plot the mean-evoked dF_F traces from all SST cells in mCherry group +/- C21
dataset = Data_ByCell(Data_ByCell.group=='mCherry',:); 

BL = []; C21 = [];

% Choice of the variable to plot HERE
for row=1:height(dataset)
    if dataset.condition(row) == "BL" 
        BL = [BL,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    elseif dataset.condition(row) == "C21" 
        C21 = [C21,dataset.allstimsallROIs_dF_mean(row)]; % change the variable
    end    
end
clear row

mean_BL = mean(cell2mat(BL'), "omitnan"); std_BL = std(cell2mat(BL'), "omitnan")/sqrt(width(BL)); % computation of the std error deviation
mean_C21 = mean(cell2mat(C21'), "omitnan"); std_C21 = std(cell2mat(C21'), "omitnan")/sqrt(width(C21));

xLabels = {'0', '1', '2', '3', '4'};
x=[0:41]*(1/10.0143333333333); %X-axis in sec (1:42 frames * 1/framerate)
means = [mean_BL', mean_C21']; stdErrors = [std_BL', std_C21'];

% Plot the bar chart with error bars for the standard errors
figure_meanEvokedResp = figure;
axes_meanEvokedResp = axes('Parent',figure_meanEvokedResp);
hold(axes_meanEvokedResp,'on');

style={'-k','-r'}; %Define colors for each whisker stimulated
greybar=NaN(size(mean_BL,2),1);
greybar(11:20)=20;
area(x,greybar,-5,'FaceColor',[.9,.9,.9],'EdgeColor','none');

for i=1:size(means,2) %error fill plot with shaded s.e.m.
    errorfill(x,means(:,i)', stdErrors(:,i)', style{i}); % parameters to modify for the design (if needed)
    hold on
end
hold off;

% Add the x-axis labels and title
set(axes_meanEvokedResp, 'XTick', [0:4], 'XTickLabel', xLabels);
xlabel('Time (sec)'); 
ylabel('Z-score'); % Choice of the name of the variable plot HERE
title('PVcre-mCherry - Mean whisker-evoked dF_F response of all cells SSTs'); % Choice of the name of the plot HERE
ylim([0 0.14]) % add a limit here if needed to scale different plots together reliability
xlim([-0.1 (max(x)+0.1)])

%Get statistics
dataset = Data_ByCell(Data_ByCell.group=='mCherry',:); 
dataset.condition = reordercats(dataset.condition, {'BL', 'C21'}); 
normality = lillietest(dataset.MI_dF) %1, not normal
lme = fitlme(dataset,'MI_dF ~ condition + (1 | MouseID)','FitMethod','REML')

clear BL C21 mean_BL mean_C21 std_BL std_C21 means stdErrors x xLabels...
    greybar style dataset i normality lme

%% Section 8 (Required): Load the spontaneous data from PV DREADD Experiments - Analysis of SST Cells
clear %clear workspace variables

load("/Users/wzeiger/Library/CloudStorage/Box-Box/ZeigerLab Sync/Lab Shared/Data/Baruc/SST_DREADDs_Cre-Flip/Spontaneous Analysis/SST_DREADDs_SpontDataByCell_dF.mat")

% Note if one wants to match up cells between spontaneous and stimulus evoked:
% Spontaneous and Stm_C1 cell selection is not aligned for m593.
% Spontaneous movies match up with the cell selection as indicated on the “m593_C1-1.tif” and “m593_C1-2.tif”.
% Stm_C1 cell selection is matched between baseline and +C21. Cells 1, 2, 7, 9, 5, 6, 3, 4 match up with spontaneous. 
%Cell 8 from spontaneous is not selected. 3 additional cells at 353,80; 233,113; and 441,356 are selected.

% Calculate data for dF_F trace
framerate=10.0143333333333;
Data=Spont_Data_ByCell;

for ii=1:height(Data)
    %Calculate the baseline of dF_F
    F=double(Data.F{ii});
    Fsize=size(F');
    baselineframes=round(10*framerate);

    %-----------Calculate Baseline F------------------------
    basecalc=zeros(1,Fsize(1)-baselineframes+1);
    for i=1:Fsize(1)-baselineframes+1 %find 10 seconds of baseline. RECALC based on framerate
        basecalc(i)=std(F(i:i+baselineframes-1)); %checks STD for all values
    end
    [~,mindex]=min(basecalc); %calculate the position of lowest deviation for baseline
    dF_std=std(Data.dF_F{ii}(mindex:mindex+baselineframes-1));
    Data.dF_std(ii)=dF_std;
    %----------------End Calculate Baseline F--------------------------------
    clear Fraw Fneu F Fsize baselineframes basecalc mindex
    
    %Calculate the mean dF/F for the entire trace
    Data.mean_dF(ii)=mean(Data.dF_F{ii});
    
    %Calculate the AUC of the entire dF/F trace
    Data.AUC_dF(ii)=trapz(Data.dF_F{ii})/framerate;

    %Peak finding for entire trace
    warning('off','signal:findpeaks:largeMinPeakHeight')
    
    %Set MinPeakHeight to 4x STD of baseline and MinPeakProminence to 3x
    %STD to match the Z-scored analysis settings
    [pks,locs,w,p]=findpeaks(Data.dF_F{ii},'MinPeakHeight',(4*dF_std),'MinPeakDistance',round(framerate*0.5),'MinPeakProminence',(3*dF_std));
    Data.Peaks_Amp_dF{ii}=pks;
    Data.Peaks_Loc_dF{ii}=locs;
    Data.Peaks_Width_dF{ii}=w;
    Data.Peaks_Prom_dF{ii}=p;

    traceduration=size(Data.dF_F{ii},2)/framerate;
    
    Data.PeaksAmp_mean_dF(ii)=mean(pks); %mean peak amplitude Z-score
    Data.PeakFreq_dF(ii)=size(locs,2)/traceduration; %peak frequency as peaks/sec
    clear pks locs w p traceduration dF_std
end

warning('on','signal:findpeaks:largeMinPeakHeight')
Spont_Data_ByCell=Data;
clear i ii

% Pool all cells by group
T1=Spont_Data_ByCell(:,{'condition','group'});
[G,mean_allcells]=findgroups(T1);
C=groupcounts(T1,{'condition','group'});
omean=@(x) mean(x,'omitnan');
omean2=@(x) mean(x,1,'omitnan');
sem=@(x) (std(x,'omitnan')/sqrt(size(x,1)));
sem2=@(x) (std(x,1,'omitnan')/sqrt(size(x,1)));

mean_allcells.TotalCells=C.GroupCount;
mean_allcells.AUC=splitapply(omean,Spont_Data_ByCell.AUC,G);
mean_allcells.AUC_SEM=splitapply(sem,Spont_Data_ByCell.AUC,G);
mean_allcells.mean_Z=splitapply(omean,Spont_Data_ByCell.mean_Z,G);
mean_allcells.mean_Z_SEM=splitapply(sem,Spont_Data_ByCell.mean_Z,G);
mean_allcells.PeaksAmp_mean=splitapply(omean,Spont_Data_ByCell.PeaksAmp_mean,G);
mean_allcells.PeaksAmp_mean_SEM=splitapply(sem,Spont_Data_ByCell.PeaksAmp_mean,G);
mean_allcells.PeakFreq=splitapply(omean,Spont_Data_ByCell.PeakFreq,G);
mean_allcells.PeakFreq_SEM=splitapply(sem,Spont_Data_ByCell.PeakFreq,G);

% For dF_F
mean_allcells.AUC_dF=splitapply(omean,Spont_Data_ByCell.AUC_dF,G);
mean_allcells.AUC_dF_SEM=splitapply(sem,Spont_Data_ByCell.AUC_dF,G);
mean_allcells.mean_dF=splitapply(omean,Spont_Data_ByCell.mean_dF,G);
mean_allcells.mean_dF_SEM=splitapply(sem,Spont_Data_ByCell.mean_dF,G);
mean_allcells.PeaksAmp_mean_dF=splitapply(omean,Spont_Data_ByCell.PeaksAmp_mean_dF,G);
mean_allcells.PeaksAmp_mean_dF_SEM=splitapply(sem,Spont_Data_ByCell.PeaksAmp_mean_dF,G);
mean_allcells.PeakFreq_dF=splitapply(omean,Spont_Data_ByCell.PeakFreq_dF,G);
mean_allcells.PeakFreq_dF_SEM=splitapply(sem,Spont_Data_ByCell.PeakFreq_dF,G);

clear T1 G C omean omean2 sem sem2 framerate
%% Fig. S13a: Plot the AUC of the spontaneous activity trace from SST cells in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= dataset.AUC_dF;
stdErrors = dataset.AUC_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_AUC_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_AUC_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('AUC (dF/F*sec)');
title('PVcre-HM3DGq - SST - AUC_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq',:);
normality = lillietest(dataset_stats.AUC_dF) %1 not normal
lme = fitlme(dataset_stats,'AUC_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats lme k
%% Fig. S13b: Plot the Peak Frequency of the spontaneous activity trace from SST cells in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= dataset.PeakFreq_dF;
stdErrors = dataset.PeakFreq_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeakFreq_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeakFreq_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('PeakFreq (Hz)');
title('PVcre-HM3DGq - SST - PeakFreq_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq',:);
normality = lillietest(dataset_stats.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_stats,'PeakFreq_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats lme k
%% Fig. S13c: Plot the Peak Amp. of the spontaneous activity trace from SST cells in Gq group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gq',:); 
means= dataset.PeaksAmp_mean_dF;
stdErrors = dataset.PeaksAmp_mean_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeaksAmp_mean_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeaksAmp_mean_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('Peak Amp (dF/F)');
title('PVcre-HM3DGq - SST - PeaksAmp_mean_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gq',:);
normality = lillietest(dataset_stats.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_stats,'PeaksAmp_mean_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats lme k
%% Fig. S13d: Plot the AUC of the spontaneous activity trace from SST cells in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:);  
means= dataset.AUC_dF;
stdErrors = dataset.AUC_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_AUC_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_AUC_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('AUC (dF/F*sec)');
title('PVcre-HM4DGi - SST - AUC_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi',:);
normality = lillietest(dataset_stats.AUC_dF) %1 not normal
lme = fitlme(dataset_stats,'AUC_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats lme k
%% Fig. S13e: Plot the Peak Freq. of the spontaneous activity trace from SST cells in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:); 
means= dataset.PeakFreq_dF;
stdErrors = dataset.PeakFreq_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeakFreq_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeakFreq_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('PeakFreq (Hz)');
title('PVcre-HM4DGi - SST - PeakFreq_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi',:);
normality = lillietest(dataset_stats.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_stats,'PeakFreq_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats k lme
%% Fig. S13f: Plot the Peak Amp. of the spontaneous activity trace from SST cells in Gi group +/- C21
dataset = mean_allcells(mean_allcells.group=='Gi',:); 
means= dataset.PeaksAmp_mean_dF;
stdErrors = dataset.PeaksAmp_mean_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeaksAmp_mean_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeaksAmp_mean_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('Peak Amp (dF/F)');
title('PVcre-HM4DGi - SST - PeaksAmp_mean_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='Gi',:);
normality = lillietest(dataset_stats.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_stats,'PeaksAmp_mean_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats k lme
%% Fig. S13g: Plot the AUC of the spontaneous activity trace from SST cells in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= dataset.AUC_dF;
stdErrors = dataset.AUC_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_AUC_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_AUC_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.AUC_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.AUC_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('AUC (dF/F*sec)');
title('PVcre-mCherry - SST - AUC_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry',:);
normality = lillietest(dataset_stats.AUC_dF) %1 not normal
lme = fitlme(dataset_stats,'AUC_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats k lme
%% Fig. S13h: Plot the Peak Freq. of the spontaneous activity trace from SST cells in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= dataset.PeakFreq_dF;
stdErrors = dataset.PeakFreq_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeakFreq_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeakFreq_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeakFreq_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeakFreq_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('PeakFreq (Hz)');
title('PVcre-mCherry - SST - PeakFreq_dF');
% ylim([0 0.8])
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry',:);
normality = lillietest(dataset_stats.PeakFreq_dF) %1 not normal
lme = fitlme(dataset_stats,'PeakFreq_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means stdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats k lme
%% Fig. S13i: Plot the Peak Amp. of the spontaneous activity trace from SST cells in mCherry group +/- C21
dataset = mean_allcells(mean_allcells.group=='mCherry',:); 
means= dataset.PeaksAmp_mean_dF;
stdErrors = dataset.PeaksAmp_mean_dF_SEM;
conditions=unique(dataset.condition);
colors=[0.3 0.3 0.3; 1 0 0];

figure_PeaksAmp_mean_dF=bar(means,'FaceColor','flat','FaceAlpha',0.5);
hold on
for k = 1:size(means,1)
    figure_PeaksAmp_mean_dF.CData(k,:) = colors(k,:);
end

dataset_BL= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='BL',:);
dataset_C21= Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry' & Spont_Data_ByCell.condition=='C21',:);

scatter(1*(ones(size(dataset_BL,1),1)),...
    dataset_BL.PeaksAmp_mean_dF,50,'ok','MarkerFaceColor','k',...
    'XJitter','randn','XJitterWidth',0.2);
scatter(2*(ones(size(dataset_C21,1),1)),...
    dataset_C21.PeaksAmp_mean_dF,50,'or','MarkerFaceColor','r',...
    'XJitter','randn','XJitterWidth',0.2);

errorbar(means, stdErrors, 'k', 'linestyle', 'none','LineWidth',2);

hold off;

% Add the x-axis labels and title
set(gca, 'XTickLabel', conditions);
xlabel('condition');
ylabel('Peak Amp (dF/F)');
title('PVcre-mCherry - SST - PeaksAmp_mean_dF');
pbaspect([1 2.13 1])

%Get statistics
dataset_stats=Spont_Data_ByCell(Spont_Data_ByCell.group=='mCherry',:);
normality = lillietest(dataset_stats.PeaksAmp_mean_dF) %1 not normal
lme = fitlme(dataset_stats,'PeaksAmp_mean_dF ~ condition + (1 | MouseID:index)','FitMethod','REML')

clear dataset means StdErrors conditions xLabels colors dataset_BL dataset_C21...
    normality h p xdata types dataset_stats k lme