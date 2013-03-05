function f_rec = unslicing(f_sliced,sl_len,tr_area,slices)
%UNSLICING  Reconstruct full signal from uniform half-overlapping slices
%   Usage:  f_rec = unslicing(f_sliced,sl_len,tr_area,slices)
%
%   Input parameters:
%         f_sliced      : Matrix containing signal slices as columns
%         sl_len        : slice length (in samples, optional)
%         tr_area       : length of each transition area (in samples,
%                         optional, default is 'floor(sl_len/16)')
%         slices        : number of slices(columns) of f_sliced (optional)
%   Output parameters:
%         f_rec      : signal resulting from the overlap-add procedure
%
%   This function performs a windowed overlap-add procedure on the columns of
%   the input matrix 'f_sliced'. The signal is assumed to have been cut with
%   a Tukey window of plateau length 'sl_len/2-tr_area' and transition areas
%   of length 'tr_area'. The window used for the overlap-add process will be
%   a smooth dual window to such a Tukey window in the sense that unaltered
%   input created by the routine 'slicing.m' will recreate the original input
%   signal.
%

% Author: Nicki Holighaus
% Date: 04.03.13

    hopsize = sl_len/2;
    L = slices*sl_len/2;

    f_rec = zeros(L,1);
    
% Optional unslicing by using a longer Tukey window
%     
%    tr_area2 = min(hopsize-tr_area,2*tr_area);
%    tw = zeros(sl_len,1);
%    tw(1+(hopsize-tr_area)/2:(3*hopsize+tr_area)/2) = ones(hopsize+tr_area,1);
%    tw([(3*hopsize+tr_area)/2+(1:tr_area2/2),(hopsize-tr_area)/2+(-tr_area2/2+1:0)]) = hannwin(tr_area2);
%
% -------------------------------------------------

% Unslicing via dual window

    tw = zeros(sl_len,1);
    dual = @(x) (1+cos(pi*x))./(1+cos(pi*x).^2);
    tw(1+floor((hopsize+tr_area)/2):floor((3*hopsize-tr_area)/2)) = ...
        ones(hopsize-tr_area,1);
    tw([1+floor((hopsize-tr_area)/2):floor((hopsize+tr_area)/2), ...
        floor((3*hopsize-tr_area)/2)+1:floor((3*hopsize+tr_area)/2)]) = ...
        dual((-tr_area:tr_area-1)/tr_area);
    
% -------------------------------------------------
    
% Windowed overlap-add procedure    
    
    f_rec([ceil(end-.5*hopsize)+1:end,1:ceil(1.5*hopsize)]) = ...
        f_rec([ceil(end-.5*hopsize)+1:end,1:ceil(1.5*hopsize)]) + ...
        circshift(f_sliced(:,1),[floor(hopsize/2),0]).*tw;
    
    for kk=2:slices-1
        f_rec(1+ceil((kk-1.5)*hopsize):ceil((kk+0.5)*hopsize)) = ...
            f_rec(1+ceil((kk-1.5)*hopsize):ceil((kk+0.5)*hopsize)) + ...
            circshift(f_sliced(:,kk),[floor(((-1)^(kk-1))*hopsize/2),0]).*tw;
    end
    
    f_rec([1+ceil((slices-1.5)*hopsize):end,1:ceil((slices+0.5)*hopsize-L)]) = ...
        f_rec([1+ceil((slices-1.5)*hopsize):end,1:ceil((slices+0.5)*hopsize-L)]) + ...
        circshift(f_sliced(:,slices),[floor(((-1)^(slices-1))*hopsize/2),0]).*tw;