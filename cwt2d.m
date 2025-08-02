function out = cwt2d(fimg, wavname, scales, angles, varargin)

%% Managing of the input
if (nargin < 4 | nargout > 1)
  error('Argument Mismatch - Check Command Line');
end

if isnumeric(wavname)
  error('''wavname'' must be a string');
end

%% Seeing if the user want a yawtb object in output or just a
%% matrix (see explanations for 'Export' above).

[export,varargin] = getopts(varargin,'export',[],1);


%% The wavelet matlab function must exist (in lower case) in the
%% wave_defs subdir.
%% Note tha the real name can be wavename or [wavname '2d'] to
%% avoid eventuel confusion with 1D wavelets. This case is checked
%% first, before the existence of 'wavname' alone.
wavname = lower(wavname);

if (exist([wavname '2d']) >= 2) 
  wavname = [ wavname '2d'];
elseif (exist(wavname) < 2)
  error(['The wavelet ''' wavname ''' or ''' wavname ...
	 '2d'' doesn''t exist!']);
end
  

if ~exist('varargin')
  varargin = {};
end

%% Keeping varargin into output for reproducibility
if (~export)
  out.extra = varargin;
end

%% Input handling
if (~all(isnumeric(scales))) | (~all(isnumeric(angles))) 
  error('scale(s) and angle(s) must be numeric');
end

if (isempty(fimg)) fimg = 0; end

%% Choice of the normalization ('getopts' is part of yawtb: see the
%% 'utils' yawtb's subdirectory)
[NormChoice,varargin] = getopts(varargin,'norm','l2');
switch lower(NormChoice)
 case 'l2'
  norm = 1;
 case 'l1'
  norm = 0;
 case 'l0';
  norm = 2;
 otherwise %% Default: the L2 normalization.
  norm = 1;
end

%% Are we in the contrast normalization case ?
[ctr,varargin] = getopts(varargin,'contrast',[],1);
if (ctr)
  if (exist([wavname '_ctr']) >= 2) 
    ctrname = [ wavname '_ctr'];
  else
    error(['The wavelet ''' wavname ''' has no contrast' ...
		    ' implemented!']);
  end    
end

%% Determining if we are in the 'Pos' mode 
[fixpos, varargin] = getopts(varargin,'pos',[]);
defexec = '';
if (~isempty(fixpos))
  if strcmp(fixpos,'inter')
    [fig, varargin] = getopts(varargin, 'fig', gcf);
    figure(fig);
    [xsel, ysel] = ginput(1);
    xsel = max(1,min(size(fimg,2),round(xsel)));
    ysel = max(1,min(size(fimg,1),round(ysel)));
    fixpos  = [xsel, ysel];
    defexec = '$cwt(ysel,xsel)';
    fprintf('CWT computed on (x:%i,y:%i)\n',xsel,ysel);
  elseif (~isnumeric(fixpos))
    error('Unknown mode for the selection of the position');
  end
  defexec = '$cwt(fixpos(2),fixpos(1))';
end

%% Determining if we are in the exec mode
[exec,varargin] = getopts(varargin, 'exec', defexec);
is_exec = ~isempty(exec);
if (is_exec)
  exec = strrep(exec,'$cwt','tmp');
  exec = strrep(exec,'$last','[out.data{end,end}]');
  exec = strrep(exec,'$rec','out.data');
  exec = strrep(exec,'$fimg','fimg');
  sep  = find(exec == ';');
  if (~isempty(sep))
    init_exec = exec(1:sep);
    exec = exec(sep+1:end);
  else
    init_exec ='';
  end
end

%% Determining the wavelet parameters
wavopts    = yawopts(varargin,wavname);

%% Let see if the user don't want a progress bar
[nopbar,varargin] = getopts(varargin, 'NoPBar', [], 1);

%% Creation of the frequency plane ('vect' is part of yawtb: see
%% utils directory)
[Hgth,Wdth] = size(fimg);
[kx,ky] = yapuls2(Wdth, Hgth);
dkxdky = abs( (kx(1,2) - kx(1,1)) * (ky(2,1) - ky(1,1)) );

nsc = length(scales);
nang = length(angles);

if (is_exec)
  out.data  = {};
  if (~isempty(init_exec))
    eval(init_exec);
  end
else
  out.data  = zeros(Hgth,Wdth,nsc,nang);
end

isloop = (nsc*nang > 1);

if (isloop)
  if (~nopbar)
    oyap = yapbar([],nsc*nang);
  end
  
  for sc = 1:nsc,
    for ang = 1:nang,
      if (~nopbar)
	oyap = yapbar(oyap,'++');
      end
      
      csc  = scales(sc);
      cang = angles(ang);
      
      [nkx, nky] = yadiro(kx, ky, csc, cang, 'freq');
      
      %% Call of the wavelet function through 'eval'.
      mask = csc^norm * feval(wavname, nkx,nky,wavopts{:});
      
      if (is_exec)
	tmp = ifft2(fimg.* conj(mask));
	out.data{sc,ang} = eval(exec);
      else
	out.data(:,:,sc,ang) = ifft2(fimg.* conj(mask));
	out.wav_norm(sc, ang) = (sum(abs(mask(:)).^2)*dkxdky)^0.5/(2*pi);
      end
      
      if ( (ctr) & (~is_exec) )
	%% The kernel in the frequency space
	mask  = csc^norm * feval(ctrname, nkx,nky,wavopts{:});
	%% The local luminence (must be real) according to this kernel
	lumin = real(ifft2(fimg.* conj(mask)));
	%% Definition of the contrast
	out.data(:,:,sc,ang) = (out.data(:,:,sc,ang) ~= 0) .* ...
	    out.data(:,:,sc,ang) ./ lumin;
	out.wav_norm(sc, ang) = (sum(abs(mask(:)).^2)*dkxdky)^0.5/(2*pi);
      end
      
    end
  end
  
  if (~nopbar)
    oyap = yapbar(oyap,'close');
  end
  
else
  
  [nkx, nky] = yadiro(kx, ky, scales, angles, 'freq');

  %% Call of the wavelet function through 'eval'.
  mask = scales^norm * feval(wavname, nkx,nky,wavopts{:});
  
  out.data = ifft2(fimg .* conj(mask));
  out.wav_norm = (sum(abs(mask(:)).^2)*dkxdky)^0.5/(2*pi);

  if (ctr)
    %% The kernel in the frequency space
    mask  = feval(ctrname, nkx,nky,wavopts{:});
    %% The local luminence (must be real) according to this kernel
    lumin = real(scales^norm * ifft2(fimg.* conj(mask)));
    %% Definition of the contrast
    out.data = (out.data ~= 0) .* out.data ./ lumin;
  end
end
  
%% Setting the output
if (export)
  out = out.data;
else
  out.type  = mfilename;
  out.wav   = wavname;
  out.para  = [ wavopts{:} ]; 
  out.sc    = scales;
  out.ang   = angles;
  %% out.extra is already recorded (see top)
  out.pos   = fixpos;
end