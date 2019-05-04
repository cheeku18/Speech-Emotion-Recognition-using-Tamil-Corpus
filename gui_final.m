function varargout = gui_final(varargin)
% GUI_FINAL M-file for gui_final.fig
%      GUI_FINAL, by itself, creates a new GUI_FINAL or raises the existing
%      singleton*.
%
%      H = GUI_FINAL returns the handle to a new GUI_FINAL or the handle to
%      the existing singleton*.
%
%      GUI_FINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINAL.M with the given input arguments.
%
%      GUI_FINAL('Property','Value',...) creates a new GUI_FINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_final_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_final_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_final

% Last Modified by GUIDE v2.5 16-Feb-2017 14:14:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_final_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_final_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_final is made visible.
function gui_final_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_final (see VARARGIN)
% Choose default command line output for gui_final
handles.output = hObject;
a=ones([200 450]);
axes(handles.axes1);imshow(a);
axes(handles.axes2);imshow(a);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_final wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_final_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in inp_voice.
function inp_voice_Callback(hObject, eventdata, handles)
% hObject    handle to inp_voice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd Databases
file=uigetfile('*.wav');
inp=wavread(file);
[ speech, fs, nbits ] = wavread(file);
cd ..
wavplay(inp,44200);
axes(handles.axes1);
plot(speech);title('Input Voice Signal');

handles.speech=speech;
handles.fs=fs;
handles.nbits=nbits;
handles.file=file;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pre_process.
function pre_process_Callback(hObject, eventdata, handles)
% hObject    handle to pre_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

speech=handles.speech;
filt_sig=medfilt2(speech,[3 3]);
wavplay(filt_sig,44200);
axes(handles.axes2);
plot(filt_sig);title('Filtered Signal');
handles.filt_sig=filt_sig;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in FFT.
function FFT_Callback(hObject, eventdata, handles)
% hObject    handle to FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filt_sig=handles.filt_sig;
[rows cols] = size(filt_sig);
fft_sig =fft(filt_sig,[rows cols]);
figure;
plot(fft_sig);title('FFT Signal');
handles.fft_sig=fft_sig;


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in features.
function features_Callback(hObject, eventdata, handles)
% hObject    handle to features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inp=handles.speech;
fs=handles.fs;
nbits=handles.nbits;
fft_sig=handles.fft_sig;

%% FFT Feature Extraction 

    f1=max(max(fft_sig));
    f2=min(min(fft_sig));
    f3=mean(mean(fft_sig));
    f5=mean(mean(abs(medfilt1(fft_sig))));
    f6=std2(fft_sig);
    p= hist(inp);
    f7= -sum(sum(p.*log2(p)));
    f8=entropy(inp,256);
   [f9, t] = FeatureTimeZeroCrossingRate(inp, 42100, 256,256);
   f9=mean(f9);
   f10=sum(sum(fft_sig));


    disp('Input FFT Features:');
    disp('FFT:Max Signal Level:');disp(f1);
    disp('FFT:Min Signal Level:');disp(f2);
    disp('FFT:Average Signal Level:');disp(f3);
    %disp('DWT:Peak Level:');disp(f4);
    disp('FFT:Median Filter Signal Level:');disp(f5);
    disp('FFT:Standard Deviation:');disp(f6);
    disp('FFT:Histogram :');disp(f7);
    disp('FFT:Entropy Level:');disp(f8);
    disp('FFT:Zero Crossing Rate:');disp(f9);
    disp('FFT:Fundamental Energy Level:');disp(f10);

% 
msgbox('FFT Features Extraction Completed');

%% MFCC Feature Extraction 

    Tw = 25;                % analysis frame duration (ms)
    Ts = 10;                % analysis frame shift (ms)
    alpha = 0.97;           % preemphasis coefficient
    M = 20;                 % number of filterbank channels 
    C = 13;                 % number of cepstral coefficients
    L = 22;                 % cepstral sine lifter parameter
   R = [ 300 3700 ];  % frequency range to consider

   
  inp = inp(:,1);

 [ MFCCs, FBEs, frames ] = ...
                    mfcc(inp, fs, Tw, Ts, alpha,@hamming, R, M, C, L );


    % Generate data needed for plotting 
    [ Nw, NF ] = size( frames );                % frame length and number of frames
    time_frames = [0:NF-1]*Ts*0.001+0.5*Nw/fs;  % time vector (s) for frames 
    time = [ 0:length(fft_sig)-1 ]/fs;           % time vector (s) for signal samples 
    logFBEs = 20*log10( FBEs );                 % compute log FBEs for plotting
    logFBEs_floor = max(logFBEs(:))-50;         % get logFBE floor 50 dB below max
    logFBEs( logFBEs<logFBEs_floor ) = logFBEs_floor; % limit logFBE dynamic range


    % Generate plots
    figure('Position', [30 30 600 800], 'PaperPositionMode', 'auto', ... 
              'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 

    subplot( 311 );
    plot( time,inp, 'k' );
    xlim( [ min(time_frames) max(time_frames) ] );
    xlabel( 'Time (s)' ); 
    ylabel( 'Amplitude' ); 
    title( 'Speech waveform'); 



subplot( 312 );
    imagesc( time_frames, [1:M], logFBEs ); 
    axis( 'xy' );
    xlim( [ min(time_frames) max(time_frames) ] );
    xlabel( 'Time (s)' ); 
    ylabel( 'Channel index' ); 
    title( 'Log (mel) filterbank energies'); 

    subplot( 313 );
    imagesc( time_frames, [1:C], MFCCs(2:end,:) ); % HTK's TARGETKIND: MFCC
    %imagesc( time_frames, [1:C+1], MFCCs );       % HTK's TARGETKIND: MFCC_0
    axis( 'xy' );
    xlim( [ min(time_frames) max(time_frames) ] );
    xlabel( 'Time (s)' ); 
    ylabel( 'Cepstrum index' );
    title( 'Mel frequency cepstrum' );

    % Set color map to grayscale
    colormap( 1-colormap('gray') ); 

  


  f11=max(max(MFCCs));
  f12=min(min(MFCCs));
  f13=mean(mean(MFCCs));
  f14=mean(mean(abs(medfilt1(MFCCs))));
  f15=std2(MFCCs);

  disp('Input MFCC Features:');
  disp('MFCC:Max Signal Level:');disp(f11);
  disp('MFCC:Min Signal Level:');disp(f12);
  disp('MFCC:Average Signal Level:');disp(f13);
  disp('MFCC:Median Filter Signal Level:');disp(f14);
  disp('MFCC:Standard Deviation:');disp(f15);   
    
    
    
qfeat=(f1+f2+f3+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15);
disp('Input Voice Signal Features:');
disp(abs(qfeat));
    
    
    
    
msgbox('Features Extraction Completed');
handles.qfeat=qfeat;
% Update handles structure
guidata(hObject, handles);







