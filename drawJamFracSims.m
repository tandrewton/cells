%% Draw sims for jamFracture.cpp
% program requires output files generated by jamFracture.cpp
% position file of deformable cells, (energy file, shape file tbd)
% output is a movie made from stitching the position file frames together
% makePlots, makeAMovie, onCluster are booleans (to be turned into flags)
%pwd should give ~/Documents/YalePhd/projects
%on cluster, pwd shoudl give /home/at965
clear;
close all;
clc;

addpath("/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/cells/bash/seq/")
%CHANGE THESE PARAMETERS AS NEEDED
N = "24";
NV = "24";
calA = "1.08";
kl = "1.0";
kb = "0";
min_seed = 1;
max_seed = 1;
att="0.5";
B=0.05;
% make plots (boolean)
makePlots = 0;
% make a movie (boolean)
makeAMovie = 1;
onCluster = 1;
dt = 0.005;

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/";

%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/fracture/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/fracture/";

if onCluster == 1
    %pwd should read /home/at965/cells
    pc_dir = "/home/at965/cells"
    subdir_pipeline = "/gpfs/loomis/project/fas/ohern/at965/cells/fracture"
    subdir_output = pc_dir+"/matlab_out"
end

mkdir(subdir_pipeline);
mkdir(subdir_output);

for att = ["0.0,0.01,0.5"]
    txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+'B='+B;
    %figure(13), clf, hold on, box on;
    for seed = min_seed:max_seed
        run_name ="fracture_N"+N+"_NV"+NV+"_calA"+calA+"_kl"+kl+"_kb"+kb+"_att"+att+"_B"+B;
        fileheader=run_name + "_seed" + seed;

        pipeline_dir =  subdir_pipeline + run_name + "/";

        output_dir = subdir_output + run_name + "/";
        mkdir(output_dir)
        
        nvestr = pipeline_dir+fileheader+'.pos.jam';
%         energystr = pipeline_dir+fileheader+'.energy';
%         energy_file = load(energystr);
%         energy_file(:,1);
%         t = energy_file(:,1);
%         K = energy_file(:,2);
%         U = energy_file(:,3);
%         E = energy_file(:,4);
%         figure(13);
%         plot(t*dt, K,'r-','linewidth',2, 'DisplayName', 'K');
%         plot(t*dt,U,'b-','linewidth',2, 'DisplayName', 'U');
%         plot(t*dt,E,'k','linewidth',2, 'DisplayName', 'K+U');
%         xlabel('$\tau$','Interpreter','latex');
%         ylabel('Energy','Interpreter','latex');
%         ax = gca;
%         ax.FontSize = 24;
%         if seed == max_seed 
%             annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%             saveas(gcf, output_dir + 'Energy'+min_seed+'_'+max_seed, 'epsc')
%         end


        % read in position data
        [trajectoryData, cell_count] = readFractureSim(nvestr);

        % get number of frames
        NFRAMES = trajectoryData.NFRAMES;
        NCELLS = trajectoryData.NCELLS;
        nv = trajectoryData.nv(1,:);
        NVTOT = sum(nv);
        L = trajectoryData.L(1,:);
        Lx = L(1);
        Ly = L(2);


        % show vertices or not
        showverts = 1;

        % get cell colors
        [nvUQ, ~, IC] = unique(nv);
        NUQ = length(nvUQ);
        cellCLR = jet(NUQ);

        % get frames to plot
        if showverts == 0
            FSTART = 1;
            FSTEP = 1;
            FEND = NFRAMES;
        else
            FSTART = 1;
            FSTEP = 1;
            FEND = NFRAMES;
        end

        if makeAMovie == 1
            moviestr = output_dir + 'frac'+fileheader+'seed_'+seed+'.mp4';
            vobj = VideoWriter(moviestr,'MPEG-4');
            vobj.FrameRate = 15;
            open(vobj);
        %end

            fnum = 1;
            figure(fnum), clf, hold on, box on;
            for ff = FSTART:FSTEP:FEND
                NCELLS = cell_count(ff);
                figure(fnum), clf, hold on, box on;
                fprintf('printing frame ff = %d/%d\n',ff,FEND);

                % get cell positions
                xpos = trajectoryData.xpos(ff,:);
                ypos = trajectoryData.ypos(ff,:);
                l0 = trajectoryData.l0(ff,:);
                for nn = 1:NCELLS
                    xtmp = xpos{nn};
                    ytmp = ypos{nn};
                    l0tmp = l0(nn);
                    clr = cellCLR(IC(nn),:);
                    if showverts == 1
                        for vv = 1:nv(nn)
                            xplot = xtmp(vv) - 0.5*l0tmp;
                            yplot = ytmp(vv) - 0.5*l0tmp;
                            for xx = -1:1
                                for yy = -1:1
                                    rectangle('Position',[xplot + xx*Lx, yplot + yy*Ly, l0tmp, l0tmp],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                                end
                            end
                        end
                    else
                        cx = mean(xtmp);
                        cy = mean(ytmp);
                        rx = xtmp - cx;
                        ry = ytmp - cy;
                        rads = sqrt(rx.^2 + ry.^2);
                        xtmp = xtmp + 0.4*l0tmp*(rx./rads);
                        ytmp = ytmp + 0.4*l0tmp*(ry./rads);
                        for xx = -1:1
                            for yy = -1:1
                                vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                                finfo = [1:nv(nn) 1];
                                patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',2);
                            end
                        end
                    end
                end

                % plot box
                plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
                axis equal;
                ax = gca;
                ax.XTick = [];
                ax.YTick = [];
                ax.XLim = [-0.25 1.25]*Lx;
                ax.YLim = [-0.25 1.25]*Ly;

                % if making a movie, save frame
                if makeAMovie == 1
                    currframe = getframe(gcf);
                    writeVideo(vobj,currframe);
                end
            end


            % close video object
            if makeAMovie == 1
                close(vobj);
            end
        end
    end
end