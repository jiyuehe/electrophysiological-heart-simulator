function varargout = assign_vertex_flag(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @assign_vertex_flag_OpeningFcn, ...
                   'gui_OutputFcn',  @assign_vertex_flag_OutputFcn, ...
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

% --- Executes just before assign_vertex_flag is made visible.
function assign_vertex_flag_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
handles.dcm_obj = datacursormode(handles.figure1);
set(handles.dcm_obj,'DisplayStyle','datatip');
movegui('center');

handles.directory = varargin{1};
handles.file_name = varargin{2};
handles.geometry = varargin{3};
handles.az = varargin{4};
handles.el = varargin{5};

handles.vertex = handles.geometry.vertex;
handles.face = handles.geometry.face;
handles.vertex_based_vertices = handles.geometry.vertex_based_vertices;
handles.edge = handles.geometry.edge;

% if exist vertex_flag then load it, if not then create it
cd(handles.directory.data_dir);
load(handles.file_name);
handles.vertex_flag = vertex_flag;
cd(handles.directory.home_dir);

handles.vertex_flag_number = 1;

handles = update_display(handles);
view(handles.axes1,handles.az,handles.el);

guidata(hObject, handles); % Update handles structure

% UIWAIT makes assign_vertex_flag wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function varargout = assign_vertex_flag_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
% delete(hObject);

function figure1_CloseRequestFcn(hObject, eventdata, handles)
% save_vertex_flag(handles);
% uiresume(handles.figure1);
delete(hObject);

function pushbutton_rotate_Callback(hObject, eventdata, handles)
rotate3d on;
function pushbutton_zoom_Callback(hObject, eventdata, handles)
zoom on;
function pushbutton_pan_Callback(hObject, eventdata, handles)
pan on;

function pushbutton_select_vertices_Callback(hObject, eventdata, handles)
set(handles.dcm_obj,'Enable','on','UpdateFcn',@myupdatefcn);
function txt = myupdatefcn(~,event_obj)
txt = {''}; % Customizes text of data tips

function pushbutton_assign_a_region_of_flags_Callback(hObject, eventdata, handles)
handles = assign_region_flags(handles);
handles = update_display(handles);
guidata(hObject, handles); % Update handles structure

function pushbutton_assign_individual_vertex_flag_Callback(hObject, eventdata, handles)
% obtain selected vertices
c_info = getCursorInfo(handles.dcm_obj);
selected_vertex_id = zeros(length(c_info),1);
for i = 1:length(c_info)
    pos = c_info(i).Position;
    
    for j = 1:size(handles.vertex,1)
        if isequal(pos,handles.vertex(j,:))
            selected_vertex_id(i) = j;
            break;
        end
    end
end
handles.vertex_flag(selected_vertex_id) = handles.vertex_flag_number;
handles = update_display(handles);
guidata(hObject, handles); % Update handles structure

function pushbutton_reset_Callback(hObject, eventdata, handles)
handles.vertex_flag = ones(size(handles.vertex,1),1);
handles = update_display(handles);
guidata(hObject, handles); % Update handles structure

function edit_vertex_flag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_vertex_flag_Callback(hObject, eventdata, handles)
handles.vertex_flag_number = str2double(get(hObject,'String'));
guidata(hObject, handles); % Update handles structure

function pushbutton_clear_pacing_sites_Callback(hObject, eventdata, handles)
vertex_flag = handles.vertex_flag;
vertex_flag(vertex_flag==2) = 1;
handles.vertex_flag = vertex_flag;
handles = update_display(handles);
guidata(hObject, handles); % Update handles structure

function pushbutton_save_Callback(hObject, eventdata, handles)
f = waitbar(0.5,'please wait...','WindowStyle','modal');
save_vertex_flag(handles);
close(f);

function pushbutton_refresh_back_face_color_Callback(hObject, eventdata, handles)
handles = update_display(handles);
guidata(hObject, handles); % Update handles structure
