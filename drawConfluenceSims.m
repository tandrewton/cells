%% Draw sims for confluence.cpp
%pwd should give ~/Documents/YalePhd/projects
clear;
close all;
clc;

%CHANGE THESE PARAMETERS AS NEEDED
N = "24";
NV = "24";
calA = "1.0";
kl = "1.0";
kb = "0";
seed = "1";
att="";
% make a movie (boolean)
makeAMovie = 0;

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/";

%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/confluence/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/confluence/";

for calA = ["1.0"]
    
    run_name ="conf_N"+N+"_NV"+NV+"_calA"+calA+"_kl"+kl+"_kb"+kb+"_att"+att;
    fileheader=run_name + "_seed" + seed;

    pipeline_dir = "/Users/AndrewTon/Documents/YalePhD/projects" + ... 
        subdir_pipeline + run_name + "/";

    output_dir = "/Users/AndrewTon/Documents/YalePhD/projects" + ... 
        subdir_output + run_name + "/";
    
    fstr = pipeline_dir+fileheader+'.pos';
    shapefstr = pipeline_dir+fileheader+'.shape';
    
    %fstr = 'Jamming/CellSim/cells/pos.test';
    %shapefstr = 'Jamming/CellSim/cells/shape.test';

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

    %plot 

    % plot fire iterations
    figure(10), clf, hold on, box on;
    plot(phi,fireit,'k-','linewidth',2);
    xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
    ylabel('FIRE its.','Interpreter','latex');
    ax = gca;
    ax.FontSize = 24;
    saveas(gcf, output_dir + 'Fire', 'epsc')

    % plot cell area over time
    figure(11), clf, hold on, box on;
    clr = winter(NCELLS);
    for nn = 1:NCELLS
        plot(phi,a(:,nn),'-','linewidth',2,'color',clr(nn,:));
    end
    xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
    ylabel('$a_\mu$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 24;
    saveas(gcf, output_dir + 'CellArea', 'epsc')

    % plot U over time
    figure(12), clf, hold on, box on;
    plot(phi,U,'bo','markersize',10);
    xlabel('$\phi_0 = L^{-2} \sum_\mu a_{0\mu}$','Interpreter','latex');
    ylabel('$U$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 24;
    saveas(gcf, output_dir + 'U', 'epsc')

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
        moviestr = output_dir + 'confluence_'+fileheader+'.mp4';
        vobj = VideoWriter(moviestr,'MPEG-4');
        vobj.FrameRate = 15;
        open(vobj);
    %end

        fnum = 1;
        figure(fnum), clf, hold on, box on;
        for ff = FSTART:FSTEP:FEND
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