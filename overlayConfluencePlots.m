%% Draw sims for confluence.cpp
%pwd should give ~/Documents/YalePhd/projects
clear;
close all;
clc;

addpath("/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/cells/bash/seq/")
%CHANGE THESE PARAMETERS AS NEEDED
N = "24";
NV = "24";
calA = "1.0";
kl = "1.0";
kb = "0";
min_seed = 1;
max_seed = 10;
att="0";
% make a movie (boolean)
makeAMovie = 0;

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/";

%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/confluence/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/confluence/";

for calA = ["1.0","1.05","1.10","1.15","1.20"]
    txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att;
    figure(9), clf, hold on, box on;
    figure(10), clf, hold on, box on;
    figure(11), clf, hold on, box on;
    figure(12), clf, hold on, box on;
    for seed = min_seed:max_seed

        run_name ="conf_N"+N+"_NV"+NV+"_calA"+calA+"_kl"+kl+"_kb"+kb+"_att"+att;
        fileheader=run_name + "_seed" + seed;

        pipeline_dir =  subdir_pipeline + run_name + "/";

        output_dir = subdir_output + run_name + "/";

        fstr = pipeline_dir+fileheader+'.pos';
        shapefstr = pipeline_dir+fileheader+'.shape';

        mkdir(output_dir)

        % read in position data
        trajectoryData = readMesophyllData(fstr);

        % get number of frames
        NFRAMES = trajectoryData.NFRAMES;
        NCELLS = trajectoryData.NCELLS;
        nv = trajectoryData.nv(1,:);
        NVTOT = sum(nv);
        L = trajectoryData.L(1,:);
        Lx = L(1);
        Ly = L(2);

        % read in shape data
        readfrmt = '%f %f %f %f %f ';
        readfrmt = [readfrmt repmat('%f %f',1,NCELLS)];
        fid = fopen(shapefstr);

        % loop over frames
        fireit = zeros(NFRAMES,1);      % number of FIRE iterations
        phi = zeros(NFRAMES,1);         % packing fraction
        P = zeros(NFRAMES,1);           % pressure
        U = zeros(NFRAMES,1);           % potential energy
        p = zeros(NFRAMES,NCELLS);      % perimeter / cells
        a = zeros(NFRAMES,NCELLS);      % area / cell
        for ff = 1:NFRAMES
            datatmp = textscan(fid,readfrmt,1);
            fireit(ff) = datatmp{2};
            phi(ff) = datatmp{3};
            P(ff) = datatmp{4};
            U(ff) = datatmp{5};
            p(ff,:) = cell2mat(datatmp(6:2:(end-1)));
            a(ff,:) = cell2mat(datatmp(7:2:end));
        end
        fclose(fid);

        % plot pressure 
        figure(9);
        plot(phi,P,'k-','linewidth',2);
        xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
        ylabel('P','Interpreter','latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed 
            annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
            saveas(gcf, output_dir + 'P'+min_seed+'_'+max_seed, 'epsc')
        end

        % plot fire iterations
        figure(10);
        plot(phi,fireit,'k-','linewidth',2);
        xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
        ylabel('FIRE its.','Interpreter','latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed 
            annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
            saveas(gcf, output_dir + 'Fire'+min_seed+'_'+max_seed, 'epsc')
        end

        % plot cell area over time
        figure(11);
        clr = winter(NCELLS);
        for nn = 1:NCELLS
            plot(phi,a(:,nn),'-','linewidth',2,'color',clr(nn,:));
        end
        xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
        ylabel('$a_\mu$','Interpreter','latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed
            annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
            saveas(gcf, output_dir + 'CellArea'+min_seed+'_'+max_seed, 'epsc')
        end

        % plot U over time
        figure(12);
        plot(phi,U,'bo','markersize',10);
        xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
        ylabel('$U$','Interpreter','latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed
            annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
            saveas(gcf, output_dir + 'U'+min_seed+'_'+max_seed, 'epsc')
        end
    end
end