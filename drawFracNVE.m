%% Draw sims for confluence.cpp
%pwd should give ~/Documents/YalePhd/projects
clear;
close all;
clc;

addpath("/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/cells/bash/seq/")
%CHANGE THESE PARAMETERS AS NEEDED
N = "12";
NV = "20";
calA = "1.08";
kl = "1.0";
kb = "0";
min_seed = 1;
max_seed = 1;
att="0.0";
% make plots (boolean)
makePlots = 0;
% make a movie (boolean)
makeAMovie = 1;

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/";

%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/fracture/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/fracture/";

for calA = ["1.08"]
    txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att;
    figure(9), clf, hold on, box on;
    figure(10), clf, hold on, box on;
    figure(11), clf, hold on, box on;
    figure(12), clf, hold on, box on;
    for seed = min_seed:max_seed

        run_name ="frac_N"+N+"_NV"+NV+"_calA"+calA+"_kl"+kl+"_kb"+kb+"_att"+att;
        fileheader=run_name + "_seed" + seed;

        pipeline_dir =  subdir_pipeline + run_name + "/";

        output_dir = subdir_output + run_name + "/";

        %fstr = pipeline_dir+fileheader+'.pos';
        %shapefstr = pipeline_dir+fileheader+'.shape';
        %fstr = 'pos.test';
        %shapefstr = 'shape.test';
        nvestr = 'pos.test.jam';

        mkdir(output_dir)

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

%         % read in shape data
%         readfrmt = '%f %f %f %f %f ';
%         readfrmt = [readfrmt repmat('%f %f',1,NCELLS)];
%         fid = fopen(shapefstr);
%         
%         % loop over frames
%         fireit = zeros(NFRAMES,1);      % number of FIRE iterations
%         phi = zeros(NFRAMES,1);         % packing fraction
%         P = zeros(NFRAMES,1);           % pressure
%         U = zeros(NFRAMES,1);           % potential energy
%         p = zeros(NFRAMES,NCELLS);      % perimeter / cells
%         a = zeros(NFRAMES,NCELLS);      % area / cell
%         for ff = 1:NFRAMES
%             datatmp = textscan(fid,readfrmt,1);
%             fireit(ff) = datatmp{2};
%             phi(ff) = datatmp{3};
%             P(ff) = datatmp{4};
%             U(ff) = datatmp{5};
%             p(ff,:) = cell2mat(datatmp(6:2:(end-1)));
%             a(ff,:) = cell2mat(datatmp(7:2:end));
%         end
%         fclose(fid);
% 
%         %plot 
%         if makePlots == 1
%             % plot pressure 
%             figure(9);
%             plot(phi,P,'k-','linewidth',2);
%             xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
%             ylabel('P','Interpreter','latex');
%             ax = gca;
%             ax.FontSize = 24;
%             if seed == max_seed 
%                 annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%                 saveas(gcf, output_dir + 'P'+min_seed+'_'+max_seed, 'epsc')
%             end
% 
%             % plot fire iterations
%             figure(10);
%             plot(phi,fireit,'k-','linewidth',2);
%             xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
%             ylabel('FIRE its.','Interpreter','latex');
%             ax = gca;
%             ax.FontSize = 24;
%             if seed == max_seed 
%                 annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%                 saveas(gcf, output_dir + 'Fire'+min_seed+'_'+max_seed, 'epsc')
%             end
% 
%             % plot cell area over time
%             figure(11);
%             clr = winter(NCELLS);
%             for nn = 1:NCELLS
%                 plot(phi,a(:,nn),'-','linewidth',2,'color',clr(nn,:));
%             end
%             xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
%             ylabel('$a_\mu$','Interpreter','latex');
%             ax = gca;
%             ax.FontSize = 24;
%             if seed == max_seed
%                 annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%                 saveas(gcf, output_dir + 'CellArea'+min_seed+'_'+max_seed, 'epsc')
%             end
% 
%             % plot U over time
%             figure(12);
%             plot(phi,U,'bo','markersize',10);
%             xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
%             ylabel('$U$','Interpreter','latex');
%             ax = gca;
%             ax.FontSize = 24;
%             if seed == max_seed
%                 annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%                 saveas(gcf, output_dir + 'U'+min_seed+'_'+max_seed, 'epsc')
%             end
%         end

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