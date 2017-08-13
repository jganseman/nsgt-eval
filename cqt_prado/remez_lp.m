function h = remez_lp(lgfiltre, bands, dev, weight)
%Programme Fortran Parks-McClellan.
% Restreint à la synthèse des passe bas
% de type 1 symétrique longueur impaire
% et type 2 symétrique longueur paire
lgrid = 16;	% Densité de grille >= 16
nbands = 2;
nodd = rem(lgfiltre,2);
if nodd==1
    nfcns=(lgfiltre+1)/2;
else
	nfcns = lgfiltre/2;
end
% Remplissage de la grille, de des et wt
[FreqGrille, taillegrilledense, des, wt]=grilledense(nfcns, bands, lgfiltre, lgrid, nodd, dev, weight);
% pondération par q(w)=cos(w/2) si longueur paire
if nodd == 0
    des = des(1:taillegrilledense)./cos(pi*FreqGrille(1:taillegrilledense));
    wt = wt(1:taillegrilledense).*cos(pi*FreqGrille(1:taillegrilledense));
end
% Extrema initiaux
iExt = extrema_init(nfcns, taillegrilledense);
% Allgrothime d'échange
[x,y,ad]=echange(nfcns, taillegrilledense, iExt, FreqGrille, des, wt);
% Calcul de la réponse impulsionnelle
h=RepImp(nfcns, ad, x, y, nodd);
return;

%______________________________________
function d = dd(k, n, m, x)
% Coefficients pour l'interpolation de Lagrange
% fonction D du Fortran
d = 1;
q = x(k);
for l=1:m
	xx = 2*(q - x(l:m:n));
    xnz=find(xx ~=0);
	d = d*prod(xx(xnz));
end
d=1/d;
return;
%_____________________________________
function [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l)	
c = ad./(cos(2*pi*FreqGrille(l))-x);  
E = (c*y'/sum(c) - des(l))*wt(l);
dtemp = nut*E - comp;
return;
%_____________________________________
function [x,y,ad]=echange(nfcns, taillegrilledense, iExt, FreqGrille, des, wt)
% Algorithme d'échange
% Le drapeau flag permmet de remplacer les goto
% du programme original. Il y a donc plusieurs copies des mêmes portions
%de code pour des conditions différentes.
comp = [];
MAXITERATIONS = 250;
deltal = -1;
nz = nfcns + 1;
nzz = nz + 1;
iter = 0;
jchnge = 1; %on initialise jchange pour le while
jet = fix((nfcns - 1)/15) + 1;
while jchnge > 0   %etiquette 100 (test 370:on boucle tant que les extrema changent)
   iter = iter + 1;
   if iter > MAXITERATIONS
      break;
   end
   x = cos(2*pi*FreqGrille(iExt(1:nz)));
   % coefficients ak de l'article
   for nn = 1:nz
       ad(nn) = dd(nn, nz, jet, x);
   end
   % Calcul de delta (dev dans le Fortran)
   add = (-1).^(0:nz-1);  % do 130 k=-k
   dnum = ad*des(iExt(1:nz))';
   dden = add*(ad./wt(iExt(1:nz)))';
   delta = dnum/dden;
   nu = 1;
   if delta > 0
      nu = -1;
   end
   delta = -nu*delta;
   y = des(iExt(1:nz)) + nu*delta*add./wt(iExt(1:nz));  %do 140 k=-k
   if delta <= deltal % OUCH
        disp('-- Non convergence --')
        disp('La synthèse du filtre est probablement incorrecte.')
      break;
   end
   deltal = delta;  %etiquette 150
   jchnge = 0;  % on se prépare à sortir si les extrema n'ont pas changés
   k1 = iExt(1);
   knz = iExt(nz);
   klow = 0;  %mettre à -1 pour programme C
   nut = -nu;
   j = 1;
   while j < nzz   % etiquette 200
      kup = iExt(j+1);
      l = iExt(j) + 1;
      nut = -nut;
      if j == 2
         y1 = comp;
      end
      comp = delta;
      flag = 1;  % test de bascule pour gérer les goto
      if l < kup   % on inverse le test goto 220
        [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
        if dtemp > 0  % on inverse le test go to 220
            comp = nut*E;
            l = l + 1;  %etiquette 210
                while l < kup  % on inverse le test go to 215
                   [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
                   if dtemp <= 0  % on inverse le test go to 215
                      break;
                   else
                      comp = nut*E;
                      l = l + 1;                       
                   end
                end    
            iExt(j) = l - 1;  % etiquette 215
            j = j + 1;
            klow = l - 1;
            jchnge = jchnge + 1; % extrema ont changé
            flag = 0; % test de bascule pour le goto 200
        end
      end
      if flag
         l = l - 2;  % etiquette 220 et 225
         while l > klow   %on inverse le test goto 250
            [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
            if dtemp > 0 || jchnge > 0 % groupe les tests goto 230 et goto 225
               break;
            end
            l = l - 1; % etiquette 235
         end
         if l <= klow  % test goto 240
            l = iExt(j) + 1; % etiquette 250
            if jchnge > 0  % test goto 215
               iExt(j) = l - 1; % on refait 215
               j = j + 1;
               klow = l - 1;
               jchnge = jchnge + 1; % extrema ont changé
            else  % sinon on a fait l-2 au lieu de l-1 (220 et 225)
               l = l + 1;
               while l < kup % test goto 260
                  [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
                  if dtemp > 0 
                         break;
                  end
                  l = l + 1; % sinon on a fait l-2 au lieu de l-1 (220 et 225)
               end
               if l < kup && dtemp > 0  % on n'est pas passé dans le while précédent
                  comp = nut*E;
                  l = l + 1;
                  while l < kup
                     [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
                     if dtemp <= 0
                        break; 
                     else
                        comp = nut*E;  
                        l = l + 1;
                     end
                  end    
                  iExt(j) = l - 1;  % on est revenu à l >= kup
	              j = j + 1;	
                  klow = l - 1;
                  jchnge = jchnge + 1; % extrema ont changé
               else  % if l < kup && dtemp > 0 n'est pas vérifié
                  klow = iExt(j);
                  j = j + 1;
               end
            end
         elseif dtemp > 0  % ici l > klow
            comp = nut*E;
            l = l - 1;
            while l > klow  % tant que l > klow on fait la boucle entre 200 et goto 200
               [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
               if dtemp <= 0
                  break;
               else
                  comp = nut*E;
                  l = l - 1;                  
               end
            end
            klow = iExt(j);
            iExt(j) = l + 1;
            j = j + 1;
            jchnge = jchnge + 1; % extrema ont changé
         else
            klow = iExt(j);
            j = j + 1;
         end
      end
   end
   while j == nzz  % etiquette 300 
      k1 = min([k1 iExt(1)]);     % if k1.GT.iExt(1)
      knz = max([knz iExt(nz)]);  % if knz.LT.iExt(nz)
      nut1 = nut;
      nut = -nu;
      l = 1;
      kup = k1;
      comp = comp*1.00001;
      luck = 1;
      flag = 1;   % test de bascule pour le goto 310
      while l < kup  % on inverse le test goto 315
         [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
         if dtemp > 0
            comp = nut*E;
            j = nzz;
            l = l + 1;
            while l < kup
               [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
               if dtemp <= 0
                  break;
               else
                  comp = nut*E;
                  l = l + 1;                  
               end
            end    
            iExt(j) = l - 1;
            j = j + 1;
            klow = l - 1;
            jchnge = jchnge + 1; % extrema ont changé
            flag = 0; % gestion bascule
            break;
         end
         l = l + 1;
      end
      if flag  % on fait 315 si bascule == 1
         luck = 6;
         l = taillegrilledense + 1; % etiquette 325
         klow = knz;
         nut = -nut1;
         comp = y1*1.00001;
         l = l - 1;
         while l > klow
            [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
            if dtemp > 0   %on inverse le test pour le break dans le if suivant
               j = nzz;
               comp = nut*E;
               luck = luck + 10;
               l = l - 1;  % etiquette 330
               while l > klow  % on inverse le test goto 340
                  [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
                  if dtemp <= 0
                     break;
                  else
                     comp = nut*E;
                     l = l - 1;
                  end
               end
               klow = iExt(j);  
               iExt(j) = l + 1;
               j = j + 1;
               jchnge = jchnge + 1; % extrema ont changé
               flag = 0; % on est passe, bascule = 0
               break;
            end
            l = l - 1;
         end
         if flag
            if luck ~= 6  % on inverse le test en 340
               iExt(1:nz) = [k1 iExt(1:nfcns)];
               jchnge = jchnge + 1; % extrema ont changé
            end
            break;
         end
      end
   end
   if j > nzz % arrive-t-on en 320 ?
      if luck > 9  % oui, et luck >9 alors goto 350
         %kn=iExt(nzz) puis iExt(nz)=kn
         iExt(1:nfcns) = iExt(2:nz);
         iExt(nz) = iExt(nzz);
         jchnge = jchnge + 1; % extrema  ont changé
      else   % oui, et luck <= 9
         y1 = max([y1 comp]);
         k1 = iExt(nzz);
         l = taillegrilledense + 1;
         klow = knz;
         nut = -nut1;
         comp = y1*1.00001;
         l = l - 1;
         while l > klow  % tant que l > klow, on fait ce qui est sous 330
            [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
            if dtemp > 0  % on inverse le test pour le break dans le if suivant
               j = nzz;
               comp = nut*E;
               luck = luck + 10;
               l = l - 1;
               while l > klow
                  [dtemp, E]= gee(FreqGrille, ad, x, y, des, wt, comp, nut, l);
                  if dtemp <= 0
                     break;
                  else
                     comp = nut*E;
                     l = l - 1;
                  end
               end
               klow = iExt(j);
               iExt(j) = l + 1;
               j = j + 1;
               jchnge = jchnge + 1; % extrema ont changé
               %kn=iExt(nzz) puis iExt(nzz)=kn
               iExt(1:nfcns) = iExt(2:nz);
               iExt(nz) = iExt(nzz);
               break;
            end
            l = l - 1;
         end
         if luck ~= 6 % on arrive en 340 
            iExt(1:nz) = [k1 iExt(1:nfcns)];
            jchnge = jchnge + 1; % extrema ont changé
         end
      end  
   end
end
return;
function  h=RepImp(nfcns, ad, x, y, nodd)
% Transformée de Fourier inverse
nz = nfcns+1;
nzz = nz+1;
nm1 = nfcns - 1;
fsh = 1.0e-6;
x(nzz) = -2;
cn = 2*nfcns - 1;
delf = 1/cn;
l = 1;
for j = 1:nfcns
   ft = (j-1)*delf;
   xt = cos(2*pi*ft);
   xe = x(l);
   % Test après 410
   % tant que xt<=xe && (xe-xt) >= fsh
   % on boucle sur l=l+1 et xe = x(l)
   while ((xt <= xe) && ((xe-xt) >= fsh))
      l = l + 1;
      xe = x(l);
   end
   if ((xt-xe) < fsh) % si (xt-xe) < fsh
      a(j) = y(l);
   else               % si (xt-xe) >= fsh
      c = ad./(cos(2*pi*ft)-x(1:nz));  
      a(j) = y*c'/sum(c);
   end
   if (l > 1) % etiquette 430
       l=l-1;
   end
end
dden = 2*pi/cn;
for j = 1:nfcns
   dnum = (j-1)*dden;
   if nm1 < 1
      alpha(j) = a(1);
   else
      alpha(j) = a(1) + 2*a(2:nfcns)*cos(dnum*(1:nm1)');
   end
end
alpha = [alpha(1) 2*alpha(2:nfcns)]/cn;
if nfcns <= 3
   alpha(nfcns + 1) = 0;
   alpha(nfcns + 2) = 0;
end

% On convertit les alpha en réponse impulsionnelle
% La réponse est symétrique

if nodd ~= 0
  h = [.5*alpha(nfcns:-1:2) alpha(1)]; % type 1 longueur impaire
else
  h = .25*[alpha(nfcns) alpha(nfcns-1:-1:2)+alpha(nfcns:-1:3)];
  h = [h 0.5*alpha(1)+0.25*alpha(2)]; % type 2 longueur paire
end
h = [h(1:nfcns-nodd) h(nfcns:-1:1)];
return;
function [FreqGrille, taillegrilledense, des, wt]=grilledense(nfcns, bands, lgfiltre, lgrid, nodd, dev, weight)
%_______________________________________________________
np=floor(.5+(nfcns*bands(2))/(.5+(bands(2)-bands(3))));
if (np==0)
    np=1;
end
% % Verifier si bands(3)=bands(4)=0.5
% % Dans ce cas ns=0 tous les max sont dans la band passante
if ((bands(3)==bands(4)) && (bands(3)==0.5))
    np=nfcns+1;
    ns=0;
else
    ns=nfcns+1-np;
    if (ns<=1)
        ns=2;
    end
end
while ((np+ns-1)*lgrid+2<=lgfiltre)
    lgrid=2*lgrid;
end
taillegrilledense=(np+ns-1)*lgrid+2;
nb=1;
if (ns==0)
    taillegrilledense=(np-1)*lgrid+2;
    delf=1/(taillegrilledense);
    kend=taillegrilledense;
else
    delf=1/(np*lgrid);
    kend=(np*lgrid)+1;
end
for k=1:kend
    FreqGrille(k)=(k-1)*bands(2)*delf;
    des(k)=dev(nb);
    wt(k)=weight(nb);
end
if (ns~=0)
    nb=2;
    delf=1/((ns-1)*lgrid);
    for k=1:((ns-1)*lgrid)+1
        FreqGrille(k+(np*lgrid)+1)=bands(3)+(k-1)*(bands(4)-bands(3))*delf;
        des(k+(np*lgrid)+1)=dev(nb);
        wt(k+(np*lgrid)+1)=weight(nb);
    end
end
% Dans les cas pair on supprime le point en 0.5 ou en 0.5-eps
if (nodd==0 && FreqGrille(taillegrilledense) > .5-delf)
    taillegrilledense = taillegrilledense - 1;
end
return;
function iExt = extrema_init(nfcns, taillegrilledense)
% iExt tableau de taille nfcns+2
temp = (taillegrilledense-1)/nfcns;
j=1:nfcns;
iExt = fix(temp*(j-1)+1);
iExt(nfcns+1) = taillegrilledense;
iExt(nfcns+2) = taillegrilledense + 1;
return;

