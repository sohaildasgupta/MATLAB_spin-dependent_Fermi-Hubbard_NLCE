% Performs the NLCE sum.
% Calls NLCE_load to read the graphs and their corresponding coefficients.
% Adds the properties together multiplying by the weight of each graph
feature('numcores')
t = 1; u=usub
order0 = 1:order_max;
mu = -u * (m - 1) / 2;
CTlist = zeros();
T_points = 500;
splits = 5;
mu_points = 500;
mu_split = 1; %from input
%Choosing mu so that density lies between 0-2 (for SU2 and 3 till
%half-filling) .. do check the densities for each N and update the mu-vals
%when necessary

mu_dict = containers.Map({2,3,4,5,6},...
    {containers.Map({2.34,7.43,10.38,33.34},{linspace(-10,0,mu_points),linspace(-20,0,mu_points),linspace(-40,0,mu_points),linspace(-80,0,mu_points)}),...
    containers.Map({2.34,7.43,10.38,33.34},{linspace(-20,0,mu_points),linspace(-40,0,mu_points),linspace(-80,0,mu_points),linspace(-160,0,mu_points)}),...
    containers.Map({2.34,7.43,10.38,33.34},{linspace(-40,0,mu_points),linspace(-80,0,mu_points),linspace(-160,0,mu_points),linspace(-320,0,mu_points)}),...
    containers.Map({2.34,7.43,10.38,33.34},{linspace(-80,0,mu_points),linspace(-160,-5,mu_points),linspace(-320,0,mu_points),linspace(-640,0,mu_points)}),...
    containers.Map({2.34,7.43,10.38,33.34},{linspace(-60,0,mu_points),linspace(-160,0,mu_points),linspace(-220,0,mu_points),linspace(-700,0,mu_points)})});

mu_points = uint16(mu_points/splits);
try
	mu_u_dict=mu_dict(m);
	mu_arr = mu_u_dict(double(u));
	if mu_split>0
		mu_array = mu_arr((mu_split-1)*mu_points+1:mu_split*mu_points);
	else
		mu_array = mu_arr;	
	end
catch
	mu_points = 1;
	mu_array = linspace(0,0,mu_points);
end

T_array = linspace(.5,10,T_points);

writematrix(mu_array,join(['NLCE_data/mu_u=' num2str(double(u)) '_N=' num2str(double(m))  '.csv'],''));
writematrix(T_array,join(['NLCE_data/T_u=' num2str(double(u)) '_N=' num2str(double(m)) '.csv'],''));
 
obs_list = {'density' 'doubleoccupancy' 'p1' 'p2' 'p3' 'p4' 'p5' 'p6' 'energy' 'entropy' 'cnn' 'ninj'};
time=tic;
    for order=order0    
        [g, coefficient] = NLCE_load(order);
        for k = 1:numel(g)
            key = key_gen(g{k});
            l = numnodes(g{k});   
            filename = join(["ED_data/",key, "t=", double(t), "u=", double(u), "mu=", double(mu), "n=", n, "m=", m, "D=", dno, 'ED.mat']);
            [side1, side2] = findedge(g{k});
            itemlist = {'density', 'doubleoccupancy','p1','p2','p3','p4','p5','p6','fileio','correlator',side1',side2'};
       
            try
                disp(['Try to find ED data for graph ', num2str(key), '...']);
                load(filename); % try to load result of the ED to one given graph
                disp(['ED data for graph ', num2str(key), ' found.']);
            catch
                disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
                tic           
                [spectra, nsigmapermute, testn, nimatrix, domatrix, p1matrix, p2matrix,p3matrix,p4matrix,p5matrix,p6matrix cnn, ninj] = ED_solver(g{k},t,u,mu,n,m,dno,itemlist{:});
                toc
                disp(['ED calculation for graph ', num2str(key), ' finished.']);
            end
        
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, "D=", dno, 'ED measurement start']));
            disp(join(["Measured item :", itemlist{:}]));   
            %Computing thermal average
            tic 
            try
				if mu_split>0
                	nlce = join(['NLCE_LDA_mat_files/',key, "t=", double(t), "u=", double(u), "mu=", double(mu), "n=", n, "m=", m, "D=", dno, 'split=',mu_split , 'LDA.mat']);
				else
                	nlce = join(['NLCE_LDA_mat_files/',key, "t=", double(t), "u=", double(u), "mu=", double(mu), "n=", n, "m=", m, "D=", dno, 'LDA.mat']);
				end
				load(nlce);
                T_check(:) = abs(alist(1:mu_points:end,1)'- T_array );
                mu_check(:) = abs(alist(1:mu_points,2)'- mu_array );
                if T_check<1e-14 
                    disp("T data match...");
                else
                    disp("T data do not match");
                    clear('alist');
                    eval(alist);
                end
                if mu_check<1e-14
                    disp("mu data match...");
                else
                    disp("mu data do not match");
                    clear('alist');
                    eval(alist);
                end
		disp("NLCE data found.");
            catch
		disp("NLCE data not found. Performing NLCE sum...");
                alist = zeros(mu_points*T_points,2+numel(obs_list));
                for index = 1:mu_points*T_points                      
                    mu_index = mod(index,mu_points) + mu_points*eq(mod(index,mu_points),0);
                    T_index = uint16(ceil(index/mu_points));
                    deltamu = mu - mu_array(mu_index);
                    T = T_array(T_index);
                    density = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, nimatrix); % particle # per site, needed to transform to per unit cell
                    doaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, domatrix);
                    p1accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p1matrix);
                    p2accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p2matrix);
                    p3accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p3matrix);
                    p4accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p4matrix);
                    p5accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p5matrix);
                    p6accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, p6matrix);
                    densityaccum = numnodes(g{k})*density; % Total number of particles in the cluster
                    energyaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, 'energy') + (mu_array(mu_index) + u * (m - 1) / 2)*densityaccum ; % entropy per site, needed to transform to per unit cell
                    entropyaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, 'entropy'); % entropy per site, needed to transform to per unit cell
                    cnnaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, cnn) ;
                    ninjaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj);                       
                    alist(index,:) = [T mu_array(mu_index) densityaccum doaccum p1accum p2accum p3accum p4accum p5accum p6accum energyaccum entropyaccum sum(cnnaccum) sum(ninjaccum)]; 
                end
                nlce = join(['NLCE_LDA_mat_files/',key, "t=", double(t), "u=", double(u), "mu=", double(mu), "n=", n, "m=", m, "D=", dno, 'LDA.mat']);
                save(nlce,'alist');
            end
            CTlist = CTlist + coefficient(k)*alist(:,3:end);
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, "D=", dno, 'ED measurement finish']));

        end
        
        path_str = join(['NLCE_data/N=' num2str(m) '/U=' num2str(u) '/order=' num2str(order)],''); 
        for i=1:numel(obs_list)
            if strcmp(obs_list{i},'cnn') || strcmp(obs_list{i},'ninj')
                if order>1
                    writematrix(reshape(CTlist(:, i),numel(mu_array),numel(T_array))', join([path_str '/' obs_list{i} '.csv'],''));        
                end
            else 
                writematrix(reshape(CTlist(:, i),numel(mu_array),numel(T_array))', join([path_str '/' obs_list{i} '.csv'],''));        
            end
        end         
    end
toc(time);
   
    
    
        
        
       
    


