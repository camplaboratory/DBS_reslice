function varargout = kjm_dbs_review(varargin)
% KJM_DBS_REVIEW MATLAB code for kjm_dbs_review.fig
%      KJM_DBS_REVIEW, by itself, creates a new KJM_DBS_REVIEW or raises the existing
%      singleton*.
%
%      H = KJM_DBS_REVIEW returns the handle to a new KJM_DBS_REVIEW or the handle to
%      the existing singleton*.
%
%      KJM_DBS_REVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KJM_DBS_REVIEW.M with the given input arguments.
%
%      KJM_DBS_REVIEW('Property','Value',...) creates a new KJM_DBS_REVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kjm_dbs_review_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kjm_dbs_review_OpeningFcn via varargin.
%
% Part of the "Hippotaxy" tool.
%     This is a program for comparing brain imaging data, and placing it into
%     Hippocampal stereotactic space, as described by the manuscript:
%     "Hippocampal stereotaxy: A novel mesial temporal stereotactic
%     coordinate system", by Kai Miller and colleagues, and is currently in
%     submission. Please cite this manuscript in any setting (manuscripts,
%     talks) where this program was used. MATLAB's GUIDE tool was used to
%     create this GUI.
%     
%     Copyright (C) 2015, Kai J Miller, Stanford Neurosurgery
%     kai.miller@stanford.edu, kjmiller@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Last Modified by GUIDE v2.5 23-Jan-2020 18:17:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kjm_dbs_review_OpeningFcn, ...
                   'gui_OutputFcn',  @kjm_dbs_review_OutputFcn, ...
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


% --- Executes just before kjm_dbs_review is made visible.
function kjm_dbs_review_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kjm_dbs_review (see VARARGIN)

% Choose default command line output for kjm_dbs_review
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kjm_dbs_review wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kjm_dbs_review_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PUT ALL OF CREATEFCN STUFF HERE %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function SagSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function CorSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function AxSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function SaveName_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function AdHocLabel_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    
    
function AxSliceEditBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    

function SagSliceEditBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    

function CorSliceEditBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end   
    
function CLimLowPre_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLimHighPre_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
function ZoomFactor_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    

function ZoomSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
function CLimLowPost_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
function CLimHighPost_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDITABLE CONTENT STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA INPUT AND OUTPUT

% --- Executes on button press in ImportMRI.
function ImportMRI_Callback(hObject, eventdata, handles)
    % start over
    clear handles.data;

    %% use spm_select to get MRI -- PRE
    [tmp]=spm_select(1,'image','select PRE mri rotated into hippo space');
    handles.data.brain_info=spm_vol(tmp); [brain_vol]=spm_read_vols(handles.data.brain_info);
    
    % scale MR intensity 0-1
    brain_vol=brain_vol/(max(reshape(brain_vol,1,[])));
    handles.clims=[0 1];
    
    % assign brain image
    handles.data.brain_vol=brain_vol;

    %% use spm_select to get MRI -- POST
    [tmp]=spm_select(1,'image','select POST image, resliced into hippo image');
    tmp2=spm_vol(tmp); [brain_vol]=spm_read_vols(tmp2);
    
    % scale MR intensity 0-1
    brain_vol=brain_vol/(max(reshape(brain_vol,1,[])));
    handles.climsPost=[0 1];
    
    % assign brain image
    handles.data.brain_volPost=brain_vol;
    
    %% set starting and max cross-sections
    handles.currentXYZ=round(size(brain_vol)/2);
    handles.MaxXYZ=size(brain_vol);

    guidata(hObject,handles);imageInitiate(hObject, eventdata, handles)

function LoadROImask_Callback(hObject, eventdata, handles) % --- Executes on button press in LoadROImask.
    %% use spm_select to get MRI -- POST
    [tmp]=spm_select(1,'image','select ROI mask -- must be same dimensions as Pre/Post MRIs');
    tmp3=spm_vol(tmp); [mask_vol]=spm_read_vols(tmp3);
    
    % scale MR intensity 0-1
    mask_vol=mask_vol>0;
    
    % check that pre/post exist
    if ~isfield(handles.data,'brain_vol'),
        error('No MRI volume','No MRI volume')
    % check that dimensions are same as pre/post 
    elseif prod(size(handles.data.brain_vol)==size(mask_vol))~=1;
        error('MRI and ROI sizes are different','MRI and ROI sizes are different')
    else
    % assign ROI volume
    handles.data.ROIVol=mask_vol;
    end
    
    guidata(hObject,handles);imageUpdate(hObject, eventdata, handles);

function SaveMRI_Callback(hObject, eventdata, handles) % --- Executes on button press in SaveMRI.
%% save data
    data=handles.data; 
    uisave({'data'})
    
function LoadMRIMarks_Callback(hObject, eventdata, handles) % --- Executes on button press in LoadMRIMarks.
%% load data
    uiopen('.mat')
    handles.data=data;
    handles.currentXYZ=round(size(handles.data.brain_vol)/2);
    handles.MaxXYZ=size(handles.data.brain_vol);
    handles.clims=[0 1];    handles.climsPost=[0 1];
    guidata(hObject,handles);imageInitiate(hObject, eventdata, handles)    

    
function ExportROInii_Callback(hObject, eventdata, handles) % --- Executes on button press in ExportROInii.
%% toggle radiobutton to save ROI volume as .nii
    if (get(hObject,'Value') == get(hObject,'Max'))
        handles.ExportROInii=1;
    else
        handles.ExportROInii=0;
    end
    guidata(hObject,handles);
    
function InvertX_Callback(hObject, eventdata, handles)% --- Executes on button press in InvertX.
%% This is a toggle radiobutton for inverting if hippocampus is on the left so that "X" is always medial-lateral (with lateral +)
    if (get(hObject,'Value') == get(hObject,'Max'))
        handles.InvertX=1;
    else
        handles.InvertX=0;
    end
    guidata(hObject,handles);

function ExportROIsPoints_Callback(hObject, eventdata, handles) % --- Executes on button press in ExportROIsPoints.
%% export data in hippocampal coordinates

    origin=round(-handles.data.brain_info.mat(1:3,4)).';
    
    %% LandMarks - recenter at origin (note, all is in reshaped 1x1x1 voxels after hippotaxy export), invert if necessary
    if sum(prod(handles.data.LMark,2)>0)>0
    ValidChans=reshape(find(prod(handles.data.LMark,2)>0),1,[]);
    disp(['LandMarks noted: ' num2str(ValidChans)])    
    LMark_export=handles.data.LMark(ValidChans,:)-ones(length(ValidChans),1)*origin;        
    if handles.InvertX==1; % Need to flip x index if it is on L (?or R) -- without having kept ac pc, will have to be manual via toggle rather than automated like in hippotaxy script
        LMark_export(:,1)=-LMark_export(:,1);
    end
    end    
    
    %% ROI volume -- get indices of points, recenter at origin, invert if necessary
    if sum(sum(sum(handles.data.ROIVol)))>0
    [x,y,z]=ind2sub(size(handles.data.ROIVol),find(handles.data.ROIVol));    
    ROI_export=[x y z]-ones(length(x),1)*origin;
    if handles.InvertX==1; % Need to flip x index if it is on L (?or R) -- without having kept ac pc, will have to be manual via toggle rather than automated like in hippotaxy script
        ROI_export(:,1)=-ROI_export(:,1);
    end
    end
    
    %% Need to deal with ad-hoc points
    if length(handles.data.ad_hoc)>0
        ad_hoc_export=handles.data.ad_hoc;
        for k=1:length(ad_hoc_export)
            ad_hoc_export(k).coordinate=ad_hoc_export(k).coordinate-origin;
            if handles.InvertX==1; % Need to flip x index if it is on L (?or R) -- without having kept ac pc, will have to be manual via toggle rather than automated like in hippotaxy script
                ad_hoc_export(k).coordinate(1)=-ad_hoc_export(k).coordinate(1);
            end
        end
    end
    %% export landmarks and ROI
    %get filename and export location from where MRI loaded, send stuff there with changed names
    MR_NamePath=handles.data.brain_info.fname;MR_NamePath([-3:0]+end)=[]; 
    %export one matlab datafile
    save([MR_NamePath '_LMarks_ROI'],'*_export');
    %export individual datafiles as delimited text files
    % L Mark
    if sum(prod(handles.data.LMark,2)>0)>0    
        fid=fopen([MR_NamePath '_LMark.txt'],'w+');
        fprintf(fid,'% 6.0f % 6.0f %6.0f\r\n',LMark_export.');
        fclose(fid);
    end
    % ROI
    if sum(sum(sum(handles.data.ROIVol)))>0   
        fid=fopen([MR_NamePath '_ROI.txt'],'w+');
        fprintf(fid,'% 6.0f % 6.0f %6.0f\r\n',ROI_export.');
        fclose(fid);
    end    
    % Ad Hoc points
    if length(handles.data.ad_hoc)>0
        for k=1:length(ad_hoc_export)
           tmp=ad_hoc_export(k).label; tmp(tmp==' ')=[];
            fid=fopen([MR_NamePath '_' tmp '.txt'],'w+');
            fprintf(fid,'% 6.0f % 6.0f %6.0f\r\n',ad_hoc_export(k).coordinate.');
            fclose(fid);
        end
    end
    
    %% export ROI volume as nii here (if desired)
    if handles.ExportROInii==1
        ROI_brain_info=handles.data.brain_info;
        ROI_brain_info.fname=[MR_NamePath '_ROIVolout.nii'];
        ROI_brain_info.pinfo(1,1)=1;
        spm_write_vol(ROI_brain_info,double(handles.data.ROIVol)); % need to scale into standard ranges
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navigating MRIs

    % Sliders
    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    % Edit Boxes - to display number of current mr slice and jump to others
    % Hints: get(hObject,'String') returns contents of CorSliceEditBox as text
    %        str2double(get(hObject,'String')) returns contents of CorSliceEditBox as a double

    
function CorSlider_Callback(hObject, eventdata, handles) % --- Executes on slider movement.
    handles.currentXYZ(2)=round(get(hObject,'Value'));
    imageUpdate(hObject, eventdata, handles);

function AxSlider_Callback(hObject, eventdata, handles)
    handles.currentXYZ(3)=round(get(hObject,'Value'));
    imageUpdate(hObject, eventdata, handles);
    
% Edit boxes to jump to designated spot
    
function CorSliceEditBox_Callback(hObject, eventdata, handles)
    tmp=round(str2double(get(hObject,'String')));
    if and(tmp>0,tmp<=handles.MaxXYZ(2)), handles.currentXYZ(2)=tmp; end
    imageUpdate(hObject, eventdata, handles);
    
function AxSliceEditBox_Callback(hObject, eventdata, handles)
    tmp=round(str2double(get(hObject,'String')));
    if and(tmp>0,tmp<=handles.MaxXYZ(3)), handles.currentXYZ(3)=tmp; end
    imageUpdate(hObject, eventdata, handles);

% --- Executes on mouse press over axes background.
function AxAxesPre_ButtonDownFcn(hObject, eventdata, handles)
    a=get(hObject); pt=round(a.Parent.CurrentPoint(1,1:2));
    handles.currentXYZ(1)=pt(1); handles.currentXYZ(2)=pt(2);
    imageUpdate(hObject, eventdata, handles);

function CorAxesPre_ButtonDownFcn(hObject, eventdata, handles)
    a=get(hObject); pt=round(a.Parent.CurrentPoint(1,1:2));
    handles.currentXYZ(1)=pt(1); handles.currentXYZ(3)=pt(2);
    imageUpdate(hObject, eventdata, handles);   

function AxAxesPost_ButtonDownFcn(hObject, eventdata, handles)
    a=get(hObject); pt=round(a.Parent.CurrentPoint(1,1:2));
    handles.currentXYZ(1)=pt(1); handles.currentXYZ(2)=pt(2);
    imageUpdate(hObject, eventdata, handles);

function CorAxesPost_ButtonDownFcn(hObject, eventdata, handles)
    a=get(hObject); pt=round(a.Parent.CurrentPoint(1,1:2));
    handles.currentXYZ(1)=pt(1); handles.currentXYZ(3)=pt(2);
    imageUpdate(hObject, eventdata, handles);       

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving points in landmarks and ad hoc

% Landmarks
function Lmark01_Callback(hObject, eventdata, handles) % --- Executes on button press in Lmark01.
    set(handles.Lmark01,'String',['LandMark 1 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(1,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLmark01_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(1,:))>0
        handles.currentXYZ=handles.data.LMark(1,:); imageUpdate(hObject, eventdata, handles);
    end

function Lmark02_Callback(hObject, eventdata, handles) % --- Executes on button press in SetAC.
    set(handles.Lmark02,'String',['LandMark 2 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(2,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLmark02_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(2,:))>0
        handles.currentXYZ=handles.data.LMark(2,:); imageUpdate(hObject, eventdata, handles);
    end
    
function Lmark03_Callback(hObject, eventdata, handles) % --- Executes on button press in SetPC.
    set(handles.Lmark03,'String',['LandMark 3 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(3,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLmark03_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(3,:))>0
        handles.currentXYZ=handles.data.LMark(3,:); imageUpdate(hObject, eventdata, handles);
    end
    
function Lmark04_Callback(hObject, eventdata, handles) % --- Executes on button press in SetML1.
    set(handles.Lmark04,'String',['LandMark 4 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(4,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLmark04_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(4,:))>0
        handles.currentXYZ=handles.data.LMark(4,:); imageUpdate(hObject, eventdata, handles);
    end

function Lmark05_Callback(hObject, eventdata, handles) % --- Executes on button press in SetML2.
    set(handles.Lmark05,'String',['LandMark 5 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(5,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLmark05_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(5,:))>0
        handles.currentXYZ=handles.data.LMark(5,:); imageUpdate(hObject, eventdata, handles);
    end
    
function Lmark06_Callback(hObject, eventdata, handles) % --- Executes on button press in Lmark1.
    set(handles.Lmark06,'String',['LandMark 6 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(6,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);
function GoToLMark06_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(6,:))>0    
        handles.currentXYZ=handles.data.LMark(6,:); imageUpdate(hObject, eventdata, handles);
    end
    
function Lmark07_Callback(hObject, eventdata, handles) % --- Executes on button press in Lmark2.
    set(handles.Lmark07,'String',['LandMark 7 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(7,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles);    
function GoToLMark07_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(7,:))>0
        handles.currentXYZ=handles.data.LMark(7,:);  imageUpdate(hObject, eventdata, handles);    
    end
    
function Lmark08_Callback(hObject, eventdata, handles) % --- Executes on button press in Lmark3.
    set(handles.Lmark08,'String',['LandMark 8 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(8,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles); 
function GoToLMark08_Callback(hObject, eventdata, handles)    
    if prod(handles.data.LMark(8,:))>0
        handles.currentXYZ=handles.data.LMark(8,:);  imageUpdate(hObject, eventdata, handles);        
    end
    
function Lmark09_Callback(hObject, eventdata, handles) % --- Executes on button press in Lmark4.
    set(handles.Lmark09,'String',['LandMark 9 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(9,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles); 
function GoToLMark09_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(9,:))>0
        handles.currentXYZ=handles.data.LMark(9,:);  imageUpdate(hObject, eventdata, handles); 
    end
    
function Lmark10_Callback(hObject, eventdata, handles)% --- Executes on button press in Lmark10.
    set(handles.Lmark10,'String',['LandMark 10 - ' int2str(handles.currentXYZ(1)) ' ' int2str(handles.currentXYZ(2)) ' ' int2str(handles.currentXYZ(3))]);
    handles.data.LMark(10,:)=handles.currentXYZ; 
    imageUpdate(hObject, eventdata, handles); 
function GoToLmark10_Callback(hObject, eventdata, handles)
    if prod(handles.data.LMark(10,:))>0
        handles.currentXYZ=handles.data.LMark(10,:);  imageUpdate(hObject, eventdata, handles); 
    end

% ad hoc point saving, with label to be determined
function AdHocPoint_Callback(hObject, eventdata, handles) % --- Executes on button press in AdHocPoint.
    if isfield(handles, 'tmp'),
        % save stuff here.
        a=length(handles.data.ad_hoc);
        handles.data.ad_hoc(a+1).label=handles.tmp;
        handles.data.ad_hoc(a+1).coordinate=handles.currentXYZ;        
    else handles.Warning='Need to define name of ad hoc save point first'; 
    end
    imageUpdate(hObject, eventdata, handles);

function AdHocLabel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of AdHocLabel as text
    handles.tmp=get(hObject,'String');
    guidata(hObject,handles);

function GoToOrigin_Callback(hObject, eventdata, handles) % --- Executes on button press in GoToOrigin.    
    handles.currentXYZ=round(-handles.data.brain_info.mat(1:3,4).');  imageUpdate(hObject, eventdata, handles); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Changing intensity range, zooming, etc

% MRI intensity range
function CLimLowPre_Callback(hObject, eventdata, handles)
    tmp=str2double(get(hObject,'String'));
    if and(tmp>=0,tmp<handles.clims(2)), handles.clims(1)=tmp; end
    imageUpdate(hObject, eventdata, handles);
    
function CLimHighPre_Callback(hObject, eventdata, handles)
    tmp=str2double(get(hObject,'String'));
    if and(tmp<=1,tmp>handles.clims(1)), handles.clims(2)=tmp; end
    imageUpdate(hObject, eventdata, handles);
    
function CLimLowPost_Callback(hObject, eventdata, handles)
    tmp=str2double(get(hObject,'String'));
    if and(tmp>=0,tmp<handles.climsPost(2)), handles.climsPost(1)=tmp; end
    imageUpdate(hObject, eventdata, handles);
    
function CLimHighPost_Callback(hObject, eventdata, handles)    
    tmp=str2double(get(hObject,'String'));
    if and(tmp<=1,tmp>handles.climsPost(1)), handles.climsPost(2)=tmp; end
    imageUpdate(hObject, eventdata, handles);

% scaling zoom factor
function ZoomFactor_Callback(hObject, eventdata, handles)
    tmp=str2double(get(hObject,'String'));
    if and(tmp>=1,tmp<6), handles.zoom=tmp; 
    else handles.Warning='Zoom Between 1-6 Only'; end
    imageUpdate(hObject, eventdata, handles);

function ZoomSlider_Callback(hObject, eventdata, handles)
    handles.zoom=get(hObject,'Value');
    imageUpdate(hObject, eventdata, handles);        
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMAGE STUFF HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageUpdate(hObject, eventdata, handles)
% this function fully updates the image 
    
    guidata(hObject,handles);
    
    % update edit boxes
    set(handles.CorSliceEditBox,'String',int2str(handles.currentXYZ(2)));
    set(handles.AxSliceEditBox, 'String',int2str(handles.currentXYZ(3)));
    
    % update warning text
    set(handles.WarningText,'String',handles.Warning); 
    handles.Warning='No Warning';
 
    guidata(hObject,handles);

    % update sliders
    set(handles.CorSlider,'Value',handles.currentXYZ(2));
    set(handles.AxSlider, 'Value',handles.currentXYZ(3));
    set(handles.ZoomSlider,'Value',handles.zoom);

    % update MRI intensity ranges
    set(handles.CLimLowPre,'String',num2str(handles.clims(1)));
    set(handles.CLimHighPre,'String',num2str(handles.clims(2)));
    
    % update MRI intensity ranges
    set(handles.CLimLowPost,'String',num2str(handles.climsPost(1)));
    set(handles.CLimHighPost,'String',num2str(handles.climsPost(2)));
    
    % update Zoom
    set(handles.ZoomFactor,'String',num2str(handles.zoom));
    
    %% determine axes for zooming
    zoomWinSize=floor(handles.MaxXYZ/handles.zoom);
    
    xtmp=handles.currentXYZ(1)+floor([-1 1]*zoomWinSize(1)/2);
    if xtmp(1)<1, xrange=[1 zoomWinSize(1)]; 
    elseif xtmp(2)>handles.MaxXYZ(1), xrange=[-zoomWinSize(1)+1 0]+handles.MaxXYZ(1); 
    else xrange=xtmp; end, clear xtmp
    
    ytmp=handles.currentXYZ(2)+floor([-1 1]*zoomWinSize(2)/2);
    if ytmp(1)<1, yrange=[1 zoomWinSize(2)]; 
    elseif ytmp(2)>handles.MaxXYZ(2), yrange=[-zoomWinSize(2)+1 0]+handles.MaxXYZ(2); 
    else yrange=ytmp; end, clear ytmp

    ztmp=handles.currentXYZ(3)+floor([-1 1]*zoomWinSize(3)/2);
    if ztmp(1)<1, zrange=[1 zoomWinSize(3)]; 
    elseif ztmp(2)>handles.MaxXYZ(3), zrange=[-zoomWinSize(3)+1 0]+handles.MaxXYZ(3); 
    else zrange=ztmp; end, clear ztmp    
    
    %% update cross-sectional images

    % y - cor - pre
    axes(handles.CorAxesPre),
    tmp=imagesc(squeeze(handles.data.brain_vol(:,handles.currentXYZ(2),:)).',handles.clims);
    hold on, plot(handles.currentXYZ(1),handles.currentXYZ(3),'r.')
    axis equal, axis off,colormap(gray), set(gca,'YDir','normal')
    set(gca,'xlim',xrange,'ylim',zrange) %zoom
    set(tmp,'ButtonDownFcn', {@CorAxesPre_ButtonDownFcn, handles},'HitTest', 'on');
    
    % y - cor - post    
    axes(handles.CorAxesPost),
    tmp=imagesc(squeeze(handles.data.brain_volPost(:,handles.currentXYZ(2),:)).',handles.climsPost);
    hold on, plot(handles.currentXYZ(1),handles.currentXYZ(3),'r.')
    axis equal, axis off,colormap(gray), set(gca,'YDir','normal')
    set(gca,'xlim',xrange,'ylim',zrange) %zoom
    set(tmp,'ButtonDownFcn', {@CorAxesPost_ButtonDownFcn, handles},'HitTest', 'on');   
    
    
    %% z - ax - pre
    axes(handles.AxAxesPre)
    tmp=imagesc(squeeze(handles.data.brain_vol(:,:,handles.currentXYZ(3))).',handles.clims);
    hold on, plot(handles.currentXYZ(1),handles.currentXYZ(2),'r.')
    axis equal, axis off,colormap(gray), set(gca,'YDir','normal')
    set(gca,'xlim',xrange,'ylim',yrange) %zoom
    set(tmp,'ButtonDownFcn', {@AxAxesPre_ButtonDownFcn, handles},'HitTest', 'on');

    % z - ax - pre
    axes(handles.AxAxesPost)
    tmp=imagesc(squeeze(handles.data.brain_volPost(:,:,handles.currentXYZ(3))).',handles.climsPost);
    hold on, plot(handles.currentXYZ(1),handles.currentXYZ(2),'r.')
    axis equal, axis off,colormap(gray), set(gca,'YDir','normal')
    set(gca,'xlim',xrange,'ylim',yrange) %zoom
    set(tmp,'ButtonDownFcn', {@AxAxesPost_ButtonDownFcn, handles},'HitTest', 'on');    
%     
% hold off


    %% add plots of ROI
    % plot cross-section of trace
    if isfield(handles,'APathCor'), % coronals
        plot(handles.CorAxesPost, handles.APathCor(:,1),handles.APathCor(:,2),'g-'), 
    end    
    
    % plot cross-section of trace
    if isfield(handles,'APathAx'), % coronals
        plot(handles.AxAxesPost, handles.APathAx(:,1),handles.APathAx(:,2),'c-'), 
    end        
    
    % Plot saved ROI on coronals
    if sum(sum(sum(handles.data.ROIVol(:,handles.currentXYZ(2),:))))>0
        B = bwboundaries(squeeze(handles.data.ROIVol(:,handles.currentXYZ(2),:)));
        for k=1:length(B)
            hold on, plot(handles.CorAxesPost,B{k}(:,1), B{k}(:,2),'y');
            hold on, plot(handles.CorAxesPre,B{k}(:,1), B{k}(:,2),'y');
        end
    end
    
    % Plot saved ROI on axials
    if sum(sum(sum(handles.data.ROIVol(:,:,handles.currentXYZ(3)))))>0
        B = bwboundaries(squeeze(handles.data.ROIVol(:,:,handles.currentXYZ(3))));
        for k=1:length(B)
            hold on, plot(handles.AxAxesPost,B{k}(:,1), B{k}(:,2),'y');
            hold on, plot(handles.AxAxesPre,B{k}(:,1), B{k}(:,2),'y');
        end
    end    
    
function imageInitiate(hObject, eventdata, handles)

    guidata(hObject,handles);
    
    % initiate sliders - images
    set(handles.CorSlider,'Min',1,'Max',handles.MaxXYZ(2),'Value',handles.currentXYZ(2),'SliderStep',[1 1]/handles.MaxXYZ(2));
    set(handles.AxSlider, 'Min',1,'Max',handles.MaxXYZ(3),'Value',handles.currentXYZ(3),'SliderStep',[1 1]/handles.MaxXYZ(3));
    
    % initiate slider - Zoom
    set(handles.ZoomSlider, 'Min',1,'Max',6,'Value',1,'SliderStep',[.1 .1]);

    % initiate total slice numbers
    set(handles.CorSliceTotal,'String',['/' num2str(handles.MaxXYZ(2)) ' Coronal'])
    set(handles.AxSliceTotal,'String',['/' num2str(handles.MaxXYZ(3)) ' Axial'])
    
    % set zoom factor to 1
    handles.zoom=1;
    
    % Warning flag stuff
    handles.Warning='No Warning';
    
    % Initiate ROI Vol if doesn't exist already
    if isfield(handles.data,'ROIVol')~=1
        handles.data.ROIVol=(0*handles.data.brain_vol)>1; %mask of ROI - logical volume
    end
    
    % Initiate ROI .nii save as off
    handles.ExportROInii=0;
    
    % Initiate invert medial/lat as off (e.g. default is right hippo or AC-PC), it is flipped if left hippo.
    handles.InvertX=0;
    
    % Initiate landmarks
    if isfield(handles.data,'LMark')~=1
    handles.data.LMark=zeros(10,3);
    else %populate landmarks with appropriate labels if they are filled
        for k=1:10
            if prod(handles.data.LMark(k,:))>0
                tmp_display=['LandMark ' num2str(k) ' - ' int2str(handles.data.LMark(k,1)) ' ' int2str(handles.data.LMark(k,2)) ' ' int2str(handles.data.LMark(k,3))];
                eval(['set(handles.Lmark' num2str(k,'%02d') ',''String'',tmp_display);'])
            end
        end        
    end
    
    % Initiate ad-hoc
    if isfield(handles.data,'ad_hoc')~=1
        handles.data.ad_hoc=[];     
    end    
    imageUpdate(hObject, eventdata, handles);

    
%% stuff for tracing out slices

% coronal
function CorPath_Callback(hObject, eventdata, handles) % --- Executes on button press in PreAxPath.
    tmp=imfreehand(handles.CorAxesPost,'Closed',1);
    handles.APathCor=getPosition(tmp);
    guidata(hObject,handles); imageUpdate(hObject, eventdata, handles);
    
function CorPathReset_Callback(hObject, eventdata, handles) % --- Executes on button press in PreAxPathReset.
    if isfield(handles,'APathCor'), % area trace
        handles=rmfield(handles,'APathCor');
    end
    guidata(hObject,handles); imageUpdate(hObject, eventdata, handles);    
    
function CorAddToVol_Callback(hObject, eventdata, handles) % --- Executes on button press in PreOpAxAddToVol.
    if isfield(handles,'APathCor'), % area trace
        tmp=poly2mask(handles.APathCor(:,2), handles.APathCor(:,1), handles.MaxXYZ(1), handles.MaxXYZ(3));
        tmp2=handles.data.ROIVol(:,handles.currentXYZ(2),:);
        tmp2(tmp)=1;
        handles.data.ROIVol(:,handles.currentXYZ(2),:)=tmp2;
    end
    guidata(hObject,handles);  imageUpdate(hObject, eventdata, handles);
    
function RmCorVolSlice_Callback(hObject, eventdata, handles) % --- Executes on button press in RmAxVolSlice.
    if isfield(handles,'APathCor'), % area trace
        tmp=poly2mask(handles.APathCor(:,2), handles.APathCor(:,1), handles.MaxXYZ(1), handles.MaxXYZ(3));
        tmp2=handles.data.ROIVol(:,handles.currentXYZ(2),:);
        tmp2(tmp)=0;
        handles.data.ROIVol(:,handles.currentXYZ(2),:)=tmp2;
    end
    guidata(hObject,handles);  imageUpdate(hObject, eventdata, handles);
    
    
% axial
function AxPath_Callback(hObject, eventdata, handles) % --- Executes on button press in AxPath.
    tmp=imfreehand(handles.AxAxesPost,'Closed',1);
    handles.APathAx=getPosition(tmp);
    guidata(hObject,handles); imageUpdate(hObject, eventdata, handles);

function AxPathReset_Callback(hObject, eventdata, handles) % --- Executes on button press in AxPathReset.
    if isfield(handles,'APathAx'), % area trace
        handles=rmfield(handles,'APathAx');
    end
    guidata(hObject,handles); imageUpdate(hObject, eventdata, handles);   

function AxAddToVol_Callback(hObject, eventdata, handles) % --- Executes on button press in AxAddToVol.
    if isfield(handles,'APathAx'), % area trace
        tmp=poly2mask(handles.APathAx(:,2), handles.APathAx(:,1), handles.MaxXYZ(1), handles.MaxXYZ(2));
        tmp2=handles.data.ROIVol(:,:,handles.currentXYZ(3));
        tmp2(tmp)=1;
        handles.data.ROIVol(:,:,handles.currentXYZ(3))=tmp2;
    end
    guidata(hObject,handles);  imageUpdate(hObject, eventdata, handles);

function RmAxVolSlice_Callback(hObject, eventdata, handles) % --- Executes on button press in RmAxVolSlice.
    if isfield(handles,'APathAx'), % area trace
        tmp=poly2mask(handles.APathAx(:,2), handles.APathAx(:,1), handles.MaxXYZ(1), handles.MaxXYZ(2));
        tmp2=handles.data.ROIVol(:,:,handles.currentXYZ(3));
        tmp2(tmp)=0;
        handles.data.ROIVol(:,:,handles.currentXYZ(3))=tmp2;
    end
    guidata(hObject,handles);  imageUpdate(hObject, eventdata, handles);  
