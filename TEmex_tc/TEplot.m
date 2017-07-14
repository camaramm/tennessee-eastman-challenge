% Plots data from TE problem simulation.  Assumes
% that simulation time is in vector "tout" and 
% plant outputs are in matrix "simout".

TEdata.xy=[tout(:) simout(:,1:41)];
TEdata.iy=cell(1,41); for i=1:41; TEdata.iy{i}=i; end
TEdata.title=cell(1,41); TEdata.ylabel=cell(1,41);
TEdata.xlabel=cell(1,41); for i=1:41; TEdata.xlabel{i}='Hours'; end
TEdata.layout='2x2';
TEdata.title{1}='A Feed';
TEdata.ylabel{1}='kscmh';
TEdata.title{2}='D Feed';
TEdata.ylabel{2}='kg/hr';
TEdata.title{3}='E Feed';
TEdata.ylabel{3}='kg/hr';
TEdata.title{4}='A and C Feed';
TEdata.ylabel{4}='kscmh';
TEdata.title{5}='Recycle Flow';
TEdata.ylabel{5}='kscmh';
TEdata.title{6}='Reactor Feed Rate';
TEdata.ylabel{6}='kscmh';
TEdata.title{7}='Reactor Pressure';
TEdata.ylabel{7}='kPa gauge';
TEdata.title{8}='Reactor Level';
TEdata.ylabel{8}='%';
TEdata.title{9}='Reactor Temperature';
TEdata.ylabel{9}='Deg C';
TEdata.title{10}='Purge Rate';
TEdata.ylabel{10}='kscmh';
TEdata.title{11}='Product Sep Temp';
TEdata.ylabel{11}='Deg C';
TEdata.title{12}='Product Sep Level';
TEdata.ylabel{12}='%';
TEdata.title{13}='Product Sep Pressure';
TEdata.ylabel{13}='kPa gauge';
TEdata.title{14}='Product Sep Underflow';
TEdata.ylabel{14}='m3/hr';
TEdata.title{15}='Stripper Level';
TEdata.ylabel{15}='%';
TEdata.title{16}='Stripper Pressure';
TEdata.ylabel{16}='kPa gauge';
TEdata.title{17}='Stripper Underflow';
TEdata.ylabel{17}='m3/hr';
TEdata.title{18}='Stripper Temp';
TEdata.ylabel{18}='Deg C';
TEdata.title{19}='Stripper Steam Flow';
TEdata.ylabel{19}='kg/h';
TEdata.title{20}='Compressor Work';
TEdata.ylabel{20}='kW';
TEdata.title{21}='Reactor Coolant Temp';
TEdata.ylabel{21}='Deg C';
TEdata.title{22}='Separator Coolant Temp';
TEdata.ylabel{22}='Deg C';
for i=23:41, TEdata.ylabel{i}='Mole %'; end
comps=['A','B','C','D','E','F','G','H'];
for i=23:28, TEdata.title{i}=['Component ',comps(i-22),' to Reactor']; end
for i=29:36, TEdata.title{i}=['Component ',comps(i-28),' in Purge']; end
for i=37:41, TEdata.title{i}=['Component ',comps(i-33),' in Product']; end

   FIg1=figure;
   plot(tout,simout(:,7));
   xlabel('Hours'); ylabel(TEdata.ylabel(7)); 
   title(TEdata.title(7));
   set(FIg1,'Units','points',...
      'CloseRequestFcn','delete([FIg, FIg1]);');
   Pos=get(gcf,'Position');
   
   % Set up GUI
   FIg=figure('Units','points',...
   'CloseRequestFcn','delete([FIg, FIg1]);',...
	'MenuBar','none',...
	'Position',[Pos(1:2)+[20 -60] 180 40],...
	'Name','Signal Selection',...
   'NumberTitle','off');			% Figure window
   CALLback=['Sig=get(gcbo,''Value''); figure(FIg1);'...
         'plot(tout,simout(:,Sig)); title(TEdata.title(Sig));',...
         'xlabel(''Hours''); ylabel(TEdata.ylabel(Sig));',...
         'figure(FIg)'];
	uicontrol('Parent',FIg,'Units','points', ...
	'callback',CALLback,...
	'Position',[20 5 150 20],'String',TEdata.title,...
	'Value',7,'Style','popup');	% Inactive block edit box

   
