function [x_eye y_eye x_hnd y_hnd touching sen1 sen2 sen3 sen4] = aux_GetCalibratedEyeHandPos(task)
global IO SETTINGS
global SETTINGS
global STATE
global IO

vd = 42;
eye_offset_x=0;eye_offset_y=0;eye_gain_x=1;eye_gain_y=1;

%% Initialize eyetracker
if SETTINGS.useVPacq && ~SETTINGS.useMouse && ~libisloaded('vpx');
    vpx_Initialize;
    if ~vpx_GetStatus(1)
        error('ViewPoint software is not running!')
    end
elseif SETTINGS.useViewAPI && ~SETTINGS.useMouse
    [pSystemInfoData, SETTINGS.pSampleData, pEventData, pAccuracyData, CalibrationData] = InitiViewXAPI();
    calllib(SETTINGS.ViewAPIlibrary, 'iV_SetLogger', int32(1), 'iViewXSDK_monkeypsych.txt')
    
    disp('Connect to iViewX (eyetracking-server)')
    ret = calllib(SETTINGS.ViewAPIlibrary, 'iV_Connect', '192.168.1.2', int32(4444), '192.168.1.1', int32(5555));
    switch ret
        case 104
            msgbox('Could not establish connection. Check if Eye Tracker is running', 'Connection Error', 'modal');
        case 105
            msgbox('Could not establish connection. Check the communication Ports', 'Connection Error', 'modal');
        case 123
            msgbox('Could not establish connection. Another Process is blocking the communication Ports', 'Connection Error', 'modal');
        case 200
            msgbox('Could not establish connection. Check if Eye Tracker is installed and running', 'Connection Error', 'modal');
        otherwise
            msgbox('Could not establish connection', 'Connection Error', 'modal');
    end
    disp('Get System Info Data')
    calllib(SETTINGS.ViewAPIlibrary, 'iV_GetSystemInfo', pSystemInfoData)
    get(pSystemInfoData, 'Value')
    
    
    % ----------------------------------------------------------------------------
    % ---- start the calibration and validation process
    % ----------------------------------------------------------------------------
    
    disp('Show Accuracy')
    calllib(SETTINGS.ViewAPIlibrary, 'iV_GetAccuracy', pAccuracyData, int32(0))
    get(pAccuracyData, 'Value')
end

%{
if SETTINGS.useParallel
    sen = double(get_sensors_state(SETTINGS.pp,SETTINGS.sensor_pins));
    sen1 = sen(1);
    sen2 = sen(2);
    sen3 = sen(3);
    sen4 = sen(4);
    % this condition is not ideal, because it's rather setup specific than
    % dependent on the number of sensor pins. Sensor pins itself is
    % misleading because it only refers to parallel port
    if numel(SETTINGS.sensor_pins) > 4 % scanner DPZ with additional pin for scanner trigger
        sen3 = ~sen3; % sen(3) = 1 if jaw is moving, sen(3) = 0 if there is no jaw motion, sen(4) always 1 because no body motion detector connected
    end
else
    sen1 = 1;
    sen2 = 1;
    sen3 = 1;
    sen4 = 1;
end
%}

%{
if SETTINGS.ai
    data = getsample(IO.ai);
else
    data=NaN;
end
%}

[x_eye y_eye] = aux_GetCalibratedEyePos(task);
[x_hnd y_hnd touching] = aux_GetCalibratedHndPos(data, task);

if SETTINGS.ai && strcmp(SETTINGS.Motion_detection_interface,'DAQ')
    sen3 = data(IO.jaw)  >= 5;
    sen4 = data(IO.body) >= 5;
end


%% get raw eye position (from 0 to 1 in ViewPoint "camera field of view")
function [x,y] = aux_GetCalibratedEyePos(task)
global SETTINGS

if SETTINGS.useVPacq
    [x,y]=vpx_GetGazePointSmoothed;
    x = eye_gain_x*x + eye_offset_x;
    y = eye_gain_y*y + eye_offset_y;
    % arcs
    x = atan(x*30/SETTINGS.vd*pi/180)*180/pi;
    y = atan(y*30/SETTINGS.vd*pi/180)*180/pi;
    % real angles (gain factor 0.5)
    %     x = atan(x/SETTINGS.vd)*180/pi;
    %     y = atan(y/SETTINGS.vd)*180/pi;
elseif SETTINGS.useViewAPI
    calllib(SETTINGS.ViewAPIlibrary, 'iV_GetSample', SETTINGS.pSampleData);
    Smp = libstruct('SampleStruct', SETTINGS.pSampleData);
    x = Smp.leftEye.gazeX*SETTINGS.screen_w_cm/SETTINGS.screen_w_pix;
    y = Smp.leftEye.gazeY*SETTINGS.screen_h_cm/SETTINGS.screen_h_pix;
    
    x = eye_gain_x*x + eye_offset_x;
    y = eye_gain_y*y + eye_offset_y;
    % arcs
    x = atan(x*30/SETTINGS.vd*pi/180)*180/pi;
    y = atan(y*30/SETTINGS.vd*pi/180)*180/pi;
elseif SETTINGS.useMouse
    %mouseposition = get(0,'PointerLocation');
    [a,b]=GetMouse(SETTINGS.window,1);
    mouseposition=[a -b];
    x_M=(mouseposition(1)-SETTINGS.screen_w_pix/2)*SETTINGS.screen_w_cm/SETTINGS.screen_w_pix;
    y_M=(mouseposition(2)+SETTINGS.screen_h_pix*(SETTINGS.screen_uh_cm)/SETTINGS.screen_h_cm)*SETTINGS.screen_h_cm/SETTINGS.screen_h_pix;
    %y_M=(mouseposition(2)-SETTINGS.screen_h_pix*(SETTINGS.screen_h_cm-SETTINGS.screen_uh_cm)/SETTINGS.screen_h_cm)*SETTINGS.screen_h_cm/SETTINGS.screen_h_pix;
    x = atan(x_M/vd)*180/pi;
    y = atan(y_M/vd)*180/pi;
else
    x = NaN;
    y = NaN;
end


%% get hand position (from 0 to 1)
function [x,y, touching] = aux_GetCalibratedHndPos(data, task)
global SETTINGS

if SETTINGS.touchscreen && SETTINGS.ai
    x = round(data(1)*SETTINGS.touchscreen_calibration.x_gain + SETTINGS.touchscreen_calibration.x_offset);
    y = round(data(2)*SETTINGS.touchscreen_calibration.y_gain + SETTINGS.touchscreen_calibration.y_offset);
    if data(1) < SETTINGS.touchscreen_calibration.x_threshold || data(2) < SETTINGS.touchscreen_calibration.y_threshold
        touching = false;
        x = nan;
        y = nan;
    else
        touching = true;
    end
    [x, y] = pix2deg_xy(x, y);
    
elseif SETTINGS.UseMouseAsTouch  % use mouse input instead of touchscreen
    [a,b, button]=GetMouse(SETTINGS.window,1);
    mouseposition=[a -b];
    x_M=(mouseposition(1)-SETTINGS.screen_w_pix/2)*SETTINGS.screen_w_cm/SETTINGS.screen_w_pix;
    y_M=(mouseposition(2)+SETTINGS.screen_h_pix*(SETTINGS.screen_uh_cm)/SETTINGS.screen_h_cm)*SETTINGS.screen_h_cm/SETTINGS.screen_h_pix;
    %y_M=(mouseposition(2)-SETTINGS.screen_h_pix*(SETTINGS.screen_h_cm-SETTINGS.screen_uh_cm)/SETTINGS.screen_h_cm)*SETTINGS.screen_h_cm/SETTINGS.screen_h_pix;
    if button(1)
        x = atan(x_M/vd)*180/pi;
        y = atan(y_M/vd)*180/pi;
        touching = button(1);
        ShowCursor(2, SETTINGS.window, 1); % 1 - crosshair, 2 - index finger, 3 - four arrows, 4 - up/down arrow, 5 - left/right arrow, 6 - sand clock, 7 - stop sign, 8 - default arrow (from Windows)
    else
        x = nan;
        y = nan;
        touching = 0;
        %HideCursor(SETTINGS.window, 1);
    end
    
else
    x = NaN;
    y = NaN;
    touching = 0;
end
