function [newM,ActP] = smotherData(Maps, per)
if  per >1;                         per  = per/100;         end 
% if isfield(Maps,'E11') == 1
if size(Maps,2) == 5
    alldata  = Maps;             clearvars Maps
    [~,Maps] = reshapeData(alldata(:,1:4));     
    x     = Maps.X1;                                 y = Maps.Y1;
    xLin  = linspace(min(Maps.X1(:)),max(Maps.X1(:)),round(length(unique(x))*per,0));
    yLin  = linspace(min(Maps.Y1(:)),max(Maps.Y1(:)),round(length(unique(y))*per,0));
    [X,Y] = meshgrid(xLin,yLin); 
    X = X(:);       Y = Y(:);
    
elseif size(Maps,2) == 9
%     Maps = sortrows(Maps,[3,1,2]);
    alldata  = Maps;             clearvars Maps
    [~,Maps ] = reshapeData(alldata(:,1:6));
    x     = Maps.X1;         y = Maps.Y1;     z = Maps.Z1;
    xLin  = linspace(min(Maps.X1(:)),max(Maps.X1(:)),round(length(unique(x))*per,0));
    yLin  = linspace(min(Maps.Y1(:)),max(Maps.Y1(:)),round(length(unique(y))*per,0));
    zLin  = linspace(min(Maps.Z1(:)),max(Maps.Z1(:)),round(length(unique(z))*per,0));
    if length(unique(z))==1
        zLin = linspace(min(Maps.Z1(:)),max(Maps.Z1(:)),length(unique(z)));
    end
    [X,Y,Z] = meshgrid(xLin,yLin,zLin); 
    X = X(:);       Y = Y(:);       Z = Z(:);
end

if     per ==1 || per ==100;    newM     = alldata;
elseif size(alldata,2) == 5
	F11 = scatteredInterpolant(x(:),y(:),alldata(:,3),'natural');
    F22 = scatteredInterpolant(x(:),y(:),alldata(:,4),'natural');
	F12 = scatteredInterpolant(x(:),y(:),alldata(:,5),'natural');
    newM(:,1) = X(:);       newM(:,2) = Y(:);
    newM(:,3) = F11(X,Y); %Evaluate FE data on grid of experimental results
	newM(:,4) = F22(X,Y); %Evaluate FE data on grid of experimental results
    newM(:,5) = F12(X,Y); %Evaluate FE data on grid of experimental results
    
elseif size(alldata,2) == 9
    newM(:,1) = X(:);       newM(:,2) = Y(:);   newM(:,3) = Z(:);
    if length(zLin)==1 %memberane 
        F11 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,4),'natural');
        F22 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,5),'natural');
        F33 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,6),'natural');
        F12 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,7),'natural');
        F13 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,8),'natural');
        F23 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,9),'natural');
        newM(:,4) = F11(X,Y); %Evaluate FE data on grid of experimental results
        newM(:,5) = F22(X,Y); %Evaluate FE data on grid of experimental results
        newM(:,6) = F33(X,Y); %Evaluate FE data on grid of experimental results
        newM(:,7) = F12(X,Y); %Evaluate FE data on grid of experimental results
        newM(:,8) = F13(X,Y); %Evaluate FE data on grid of experimental results
        newM(:,9) = F23(X,Y); %Evaluate FE data on grid of experimental results
    else
        F11 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,4),'natural');
        F22 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,5),'natural');
        F33 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,6),'natural');
        F12 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,7),'natural');
        F13 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,8),'natural');
        F23 = scatteredInterpolant(alldata(:,1),alldata(:,2),alldata(:,3),alldata(:,9),'natural');
        newM(:,4) = F11(X,Y,Z); %Evaluate FE data on grid of experimental results
        newM(:,5) = F22(X,Y,Z); %Evaluate FE data on grid of experimental results
        newM(:,6) = F33(X,Y,Z); %Evaluate FE data on grid of experimental results
        newM(:,7) = F12(X,Y,Z); %Evaluate FE data on grid of experimental results
        newM(:,8) = F13(X,Y,Z); %Evaluate FE data on grid of experimental results
        newM(:,9) = F23(X,Y,Z); %Evaluate FE data on grid of experimental results
    end 
end
    ActP = round((1-length(newM)/length(alldata))*100,0);
    warning on; 
    warning(['Data is reduced by ' num2str(ActP) '%' ]); 
end	