function varargout = TAU(varargin)
% TAU MATLAB code for TAU.fig
%      TAU, by itself, creates a new TAU or raises the existing
%      singleton*.
%
%      H = TAU returns the handle to a new TAU or the handle to
%      the existing singleton*.
%
%      TAU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TAU.M with the given input arguments.
%
%      TAU('Property','Value',...) creates a new TAU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TAU_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TAU_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TAU

% Last Modified by GUIDE v2.5 17-Feb-2016 00:28:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TAU_OpeningFcn, ...
                   'gui_OutputFcn',  @TAU_OutputFcn, ...
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


% --- Executes just before TAU is made visible.
function TAU_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TAU (see VARARGIN)

% Choose default command line output for TAU
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%addlistener(handles.slider_k,'Value','PreSet',@(~,~)disp(get(handles.slider_k,'Value')));
PreSetListener = 1; % Слушатели события передвижения ползунка в реальном времени
if PreSetListener == 1
    addlistener(handles.slider_k,'Value','PreSet',@(~,~)SliderChanges(handles));
    addlistener(handles.slider_t,'Value','PreSet',@(~,~)SliderChanges(handles));
    addlistener(handles.slider_t1,'Value','PreSet',@(~,~)SliderChanges(handles));
    addlistener(handles.slider_t2,'Value','PreSet',@(~,~)SliderChanges(handles));
    addlistener(handles.slider_ksi,'Value','PreSet',@(~,~)SliderChanges(handles));
    addlistener(handles.slider_tau,'Value','PreSet',@(~,~)SliderChanges(handles)); 
    
    addlistener(handles.slider_all_y,'Value','PreSet',@(~,~)slider_all_y_Move(handles));
    addlistener(handles.slider_all_x,'Value','PreSet',@(~,~)slider_all_x_Move(handles));
end


% UIWAIT makes TAU wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TAU_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Clear_All(handles)
cla(handles.axes1); %очищение осей
cla(handles.axes2);
cla(handles.axes3);
cla(handles.axes4);

% --- Executes on button press in pushbutton_clear_all.
function pushbutton_clear_all_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Clear_All(handles);
FindOut(handles)

% --- Executes on slider movement.
function slider_k_Callback(hObject, eventdata, handles)
% hObject    handle to slider_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles)
getresult(handles)

function SliderChanges(handles)
SetSliders(handles)
if get(handles.radiobutton5, 'Value') == 1 %При перемещении ползунка
    SetAxes(handles)
    getresult(handles); 
end

function Draw(handles)
zero_koef = 0.00001; % на что заменяется ноль или коэффициенты, если на него не делится
%Отрисовка на графиках выбранного звена 
global Name Num K T T1 T2 ksi tau w Wsis A Fi Ws Ws_p y t
contents = get(handles.popupmenu_z,'String'); 
Name = contents{get(handles.popupmenu_z, 'Value')}; % Получение названия
Num = get(handles.popupmenu_z, 'Value');
K = get(handles.slider_k, 'Value');
T = get(handles.slider_t, 'Value');
ksi = get(handles.slider_ksi, 'Value');
tau = get(handles.slider_tau, 'Value');
%Частота
w_start = str2double(get(handles.edit_w_start, 'String')); %Начало
w_end = str2double(get(handles.edit_w_end, 'String')); %Конец
w_step = str2double(get(handles.edit_w_step, 'String')); %Шаг
w = w_start:w_step:w_end; %Изменение частоты

% Далее в блоке switch разные варианты звеньев. 
% Здесь case <цифра> - это номер звена из выпадающего списка на форме начиная с 1
% Нужно задать свои Wsis, A, Fi, Ws, "Ws_p = Ws(w=0)" 
% Так же нужно не забыть в выражения включить запаздывание tau
% Лучшим примером является "Апериодическое звено первого порядка" (case 3)
% Если график не зависит от частоты, то см. пример "Усилительное звено" (case 1)
% Если на графике ФЧХ есть "шипы", то см. пример "Апериодическое звено второго порядка" (case 10)
% Если на кривой разгона должен быть импульс (дельта функция), то см. строку 243
% Если при одном из коэффициентов, почти равном 0, функция не существует, то см. пример "Реальное дифференцирующее звено" (case 8) (как он отрисовывается)

switch Num
    case 1
        % Усилительное звено
        Wsis = tf(K, 'OutputDelay', tau);   %Передаточная
%A(1:length(w)) = K; - так задаются функции, если они не зависят от w
        A(1:length(w)) = K; %АЧХ
        Fi(1:length(w)) = -w * tau; %ФЧХ
        Ws(1:length(w)) = K; %АФХ
        Ws_p(1:length(w)) = K; %АФХ начальная точка
	case 2
        % Звено чистого транспортного запаздывания
        Wsis = tf(1, 'OutputDelay', tau); %Передаточная 
        %A(1:length(w)) = K; - так задаются функции, если они не зависят от w 
        A(1:length(w)) = 1; %АЧХ 
        Fi(1:length(w)) = -w*tau; %ФЧХ 
        Ws = (1.*exp(-tau.*(1i.*w))); %АФХ 
        Ws_p = (1.*exp(-tau.*(1i.*0))); %АФХ начальная точка
    case 3
        % Апериодическое звено первого порядка
        Wsis = tf(K, [T 1], 'OutputDelay', tau);   %Передаточная
        A = K./(sqrt(T^2.*w.^2+1)); %АЧХ
        Fi = -tau.*w-atan(T^2.*w.^2);   %ФЧХ
        Ws = (K.*exp(-tau.*(1i.*w)))./(T.*(1i.*w)+1); %АФХ
        Ws_p = (K.*exp(-tau.*(1i.*0)))./(T.*(1i.*0)+1); %АФХ начальная точка
	case 4
        % Идеальное интегрирующее звено 
        if T == 0
            T = zero_koef;
        end
        Wsis = tf(K, [T 0], 'OutputDelay', tau); %Передаточная 
        A = K./(T.*w); %АЧХ 
        Fi = -tau.*w-pi/2; %ФЧХ 
        Ws = (K.*exp(-tau.*(1i.*w)))./(T.*(1i.*w)); %АФХ 
        Ws_p = (K.*exp(-tau.*(1i.*0)))./(T.*(1i.*0)); %АФХ начальная точка
	case 5
        % Реальное интегрирующее звено 
        Wsis = tf(K, [T 1 0], 'OutputDelay', tau); %Передаточная 
        A = K./(w.*sqrt(T^2.*w.^2+1)); %АЧХ 
        Fi = -tau.*w-pi/2-atan(T.*w); %ФЧХ 
        Ws = (K.*exp(-tau.*(1i.*w)))./((1i.*w).*(T.*(1i.*w)+1)); %АФХ 
        Ws_p = (K.*exp(-tau.*(1i.*0)))./((1i.*0).*(T.*(1i.*0)+1)); %АФХ начальная точка
    case 6
        % Пропорционально-интегральное звено 
        if T == 0
            T = zero_koef;
        end
        Wsis = tf([K*T K], [T 0], 'OutputDelay', tau); %Передаточная 
        A = K./(T.*w).*(sqrt(T^2.*w.^2+1)); %АЧХ 
        Fi = -tau.*w-pi/2+atan(T.*w); %ФЧХ 
        Ws = (K+K./(T.*(1i.*w))).*exp(-tau.*w); %АФХ 
        Ws_p = (K+K./(T.*(1i.*0))).*exp(-tau.*0); %АФХ начальная точка
	case 7
        % Идеальное дифференцирующее звено 
        %ПХ - дельта-функция 
        Wsis = tf([K 0], 1, 'OutputDelay', tau); 
        A = K.*w; %АЧХ 
        Fi = -tau.*w+pi/2; %ФЧХ 
        Ws = K.*(1i.*w).*exp(-tau.*w); %АФХ 
        Ws_p = K.*(1i.*0).*exp(-tau.*0); %АФХ начальная точка
    case 8
        % Реальное дифференцирующее звено (дифференцирующее инерционное звено)
        Wsis = tf([K 0],[T 1],'OutputDelay', tau);   	%Передаточная
        A = K.*w./(sqrt(T^2.*w.^2+1)); 			%АЧХ
        Fi = pi/2-atan(T.*w)-tau.*w;   			%ФЧХ
        Ws =((K.*T.*w.^2)./(T^2.*w.^2+1)+(K.*(1i.*w))./(T^2.*w.^2+1)).*(exp(-tau.*(1i.*w))); %АФХ
        Ws_p =((K.*T.*0.^2)./(T^2.*0.^2+1)+(K.*(1i.*0))./(T^2.*0.^2+1)).*(exp(-tau.*(1i.*0))); %АФХ начальная точка
    case 9
        % Пропорционально-дифференциальное звено
        Wsis = tf([K*T K], 1, 'OutputDelay', tau);   %Передаточная. 
        A = K.*(sqrt(T^2.*w.^2+1)); %АЧХ
        Fi = -tau.*w+atan(T.*w);   %ФЧХ
        Ws = (exp(-tau.*(1i.*w))).*(K.*T.*1i.*w + K); %АФХ
        Ws_p = (exp(-tau.*(1i.*0))).*(K.*T.*1i.*0 + K); %АФХ начальная точка
    case 10
        % Апериодическое звено второго порядка
        Wsis = tf(K, [T^2 2*ksi*T 1], 'OutputDelay', tau); %Передаточная
        A = K./(sqrt((1-T^2.*w.^2).^2+(2.*ksi.*T.*w).^2)); %АЧХ
        Fi = -tau.*w-atan((2.*T.*ksi.*w)./(1-T^2.*w.^2)); %ФЧХ
        for i = 1:length(w) % Этот блок for нужен если у ФЧХ творится фигня с шипами
            if (w(i) > 1/T) && (tau >= 0) && (T ~= 0) % Здесь, если нужно, надо менять только условие на w(i). Но тоже не факт
                Fi(i) = Fi(i) - pi; % Точки  ФЧХ смещаются на pi вниз
            end
        end
        Ws = (K.*exp(-tau.*(1i.*w)))./((-T^2.*w.^2)+T.*ksi.*(1i.*w)+1); %АФХ
        Ws_p = (K.*exp(-tau.*(1i.*0)))./((-T^2.*0.^2)+T.*ksi.*(1i.*0)+1); %АФХ начальная точка
	case 11
        % Интегро-дифференцирующее звено
        T1 = get(handles.slider_t1, 'Value');
        T2 = get(handles.slider_t2, 'Value');
        Wsis = tf([K*T1 K], [T2 1], 'OutputDelay', tau);
        A = K.* sqrt(T1^2.*w.^2+1)./(sqrt(T2^2.*w.^2+1));
        Fi = -tau.*w-atan(T2^2.*w.^2)+ atan(T1^2.*w.^2);
        Ws = (K.*exp(-tau.*(1i.*w)).*(T1.*(1i.*w)+1))./(T2.*(1i.*w)+1);
        Ws_p = (K.*exp(-tau.*(1i.*0)).*(T1.*(1i.*0)+1))./(T2.*(1i.*0)+1);
    otherwise
        msgbox('Звено не определено!', 'Ошибка');
        %disp('Звено не определено') % вывод сообщения в рабочую область
end

%Отрисовка кривой разгона и установка времени
imp_koef = 3; % во сколько раз больше K величина импульса на графике
switch Num % switch и разные case для функций с импульсом (или которые строятся нестандартно)
    case 7
        % Идеальное дифференцирующее звено
        step_time = str2double(get(handles.edit_step_time, 'String')); 
        t=1:1:step_time; 
        y=tau.*t./t; 
        y2=0:10:(step_time-1)*10;
    case 8
        % Реальное дифференцирующее звено (дифференцирующее инерционное звено)
        if T ~= 0
            if get(handles.checkbox_kr_st_auto, 'Value') == 1 
                [y, t] = step(Wsis); %Кривая разгона в авто
            else
                step_time = str2double(get(handles.edit_step_time, 'String'));
                [y, t] = step(Wsis, step_time); %Кривая разгона с заданным временем
            end
        else
            if get(handles.checkbox_kr_st_auto, 'Value') == 1 
                [y, t] = step(tf([K 0],[zero_koef 1],'OutputDelay', tau)); %Кривая разгона в авто
            else
                step_time = str2double(get(handles.edit_step_time, 'String'));
                [y, t] = step(tf([K 0],[zero_koef 1],'OutputDelay', tau), step_time); %Кривая разгона с заданным временем
            end
        end
    case 9
        % Пропорционально-дифференциальное звено
        % Строим степ от простого K, чтобы не портить заданную выше Wsis
        if get(handles.checkbox_kr_st_auto, 'Value') == 1 
            [y, t] = step(tf(K, 1, 'OutputDelay', tau)); %Кривая разгона в авто
        else
            step_time = str2double(get(handles.edit_step_time, 'String'));
            [y, t] = step(tf(K, 1, 'OutputDelay', tau), step_time); %Кривая разгона с заданным временем
        end
        % Затем заменяем первое значение в матрице выходной величины на 0
        if tau ~= 0 % но если есть tau, то всё сдвигается
            imp = find(y > 0, 1); % индекс импульса из расчёта первого неравного нулю значения
            y(imp, 1) = K * imp_koef; % скачок с тау. y(строка, столбец)
        else % а если tau нет, то строим прямо из нуля
            y(1, 1)= K * imp_koef; % скачок из нуля
        end
    case 11
        % Интегро-дифференцирующее звено
        if T2 ~= 0
            if get(handles.checkbox_kr_st_auto, 'Value') == 1 
                [y, t] = step(Wsis); %Кривая разгона в авто
            else
                step_time = str2double(get(handles.edit_step_time, 'String'));
                [y, t] = step(Wsis, step_time); %Кривая разгона с заданным временем
            end
        end
        if (T2 == 0) && (T1 ~= 0)
            if get(handles.checkbox_kr_st_auto, 'Value') == 1 
                [y, t] = step(tf([K*T1 K], [zero_koef 1], 'OutputDelay', tau)); %Кривая разгона в авто
            else
                step_time = str2double(get(handles.edit_step_time, 'String'));
                [y, t] = step(tf([K*T1 K], [zero_koef 1], 'OutputDelay', tau), step_time); %Кривая разгона с заданным временем
            end
        end
        if (T1 == 0) && (T2 == 0)
            if get(handles.checkbox_kr_st_auto, 'Value') == 1 
                [y, t] = step(tf(K, 1, 'OutputDelay', tau)); %Кривая разгона в авто
            else
                step_time = str2double(get(handles.edit_step_time, 'String'));
                [y, t] = step(tf(K, 1, 'OutputDelay', tau), step_time); %Кривая разгона с заданным временем
            end
            % Затем заменяем первое значение в матрице выходной величины на 0
            if tau ~= 0 % но если есть tau, то всё сдвигается
                imp = find(y > 0, 1); % индекс импульса из расчёта первого неравного нулю значения
                y(imp, 1) = K * imp_koef; % скачок с тау. y(строка, столбец)
            else % а если tau нет, то строим прямо из нуля
                y(1, 1)= K * imp_koef; % скачок из нуля
            end
        end
    otherwise % Построение переходной по умолчанию
        if get(handles.checkbox_kr_st_auto, 'Value') == 1 
            [y, t] = step(Wsis); %Кривая разгона в авто
        else
            step_time = str2double(get(handles.edit_step_time, 'String'));
            [y, t] = step(Wsis, step_time); %Кривая разгона с заданным временем
        end
end

%Вывод результатов на графики
plot(handles.axes1, w, A, 'LineWidth', 2);
plot(handles.axes2, w, Fi, 'LineWidth', 2);
plot(handles.axes3, real(Ws), imag(Ws), '-', 'LineWidth', 2);
if get(handles.checkbox_afh_ps, 'Value') == 1    %Показывать точки на АФХ
    hold(handles.axes3, 'on');
    plot(handles.axes3, real(Ws), imag(Ws), '.', 'LineWidth', 2);
end
if get(handles.checkbox_afh_p, 'Value') == 1    %Начальная точка на АФХ
    hold(handles.axes3, 'on');
    plot(handles.axes3, real(Ws_p), imag(Ws_p), 'Marker','o');
end
switch Num 
    case 7
        plot(handles.axes4, y, y2, 'LineWidth', 2);
    otherwise
        plot(handles.axes4, t, y, 'LineWidth', 2);
end

grid(handles.axes1, 'on'); %сетка
grid(handles.axes2, 'on'); %сетка
grid(handles.axes3, 'on'); %сетка
grid(handles.axes4, 'on'); %сетка

function SetSliders(handles)
%Настройка слайдеров
%n_k = round(n_k, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
SliderStep_k1 = str2double(get(handles.edit_step_k1, 'String')); % при нажатии на стрелку какими частями двигаться
SliderStep_k2 = str2double(get(handles.edit_step_k2, 'String')); % при передвижении ползунка? Толщина самого ползунка (бегунка)?
if SliderStep_k1 > 1
    set(handles.edit_step_k1, 'String', '1');
    SliderStep_k1 = 1;
end

set(handles.slider_k, 'Min', str2double(get(handles.edit_k_min, 'String')))
set(handles.slider_k, 'Max', str2double(get(handles.edit_k_max, 'String')))
set(handles.slider_k, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_k_max, 'String')) SliderStep_k2/str2double(get(handles.edit_k_max, 'String'))])
n_k = get(handles.slider_k,'Value');
set(handles.text2,'String', n_k)

%n_t = round(n_t, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
set(handles.slider_t, 'Min', str2double(get(handles.edit_t_min, 'String')))
set(handles.slider_t, 'Max', str2double(get(handles.edit_t_max, 'String')))
set(handles.slider_t, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_t_max, 'String')) SliderStep_k2/str2double(get(handles.edit_t_max, 'String'))])
n_t = get(handles.slider_t,'Value');
set(handles.text3,'String', n_t)

%n_t1 = round(n_t1, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
set(handles.slider_t1, 'Min', str2double(get(handles.edit_t1_t2_min, 'String')))
set(handles.slider_t1, 'Max', str2double(get(handles.edit_t1_t2_max, 'String')))
set(handles.slider_t1, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_t1_t2_max, 'String')) SliderStep_k2/str2double(get(handles.edit_t1_t2_max, 'String'))])
n_T1 = get(handles.slider_t1,'Value');
set(handles.text15,'String', n_T1)

%n_t1 = round(n_t1, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
set(handles.slider_t2, 'Min', str2double(get(handles.edit_t1_t2_min, 'String')))
set(handles.slider_t2, 'Max', str2double(get(handles.edit_t1_t2_max, 'String')))
set(handles.slider_t2, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_t1_t2_max, 'String')) SliderStep_k2/str2double(get(handles.edit_t1_t2_max, 'String'))])
n_T2 = get(handles.slider_t2,'Value');
set(handles.text16,'String', n_T2)

%n_ksi = round(n_ksi, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
set(handles.slider_ksi, 'Min', str2double(get(handles.edit_ksi_min, 'String')))
set(handles.slider_ksi, 'Max', str2double(get(handles.edit_ksi_max, 'String')))
set(handles.slider_ksi, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_ksi_max, 'String')) SliderStep_k2/str2double(get(handles.edit_ksi_max, 'String'))])
n_ksi = get(handles.slider_ksi,'Value');
set(handles.text4,'String', n_ksi)

%n_tau = round(n_tau, 2); %Округление. Работает на версии матлаба R2015a 8.5.0.197613
set(handles.slider_tau, 'Min', str2double(get(handles.edit_tau_min, 'String')))
set(handles.slider_tau, 'Max', str2double(get(handles.edit_tau_max, 'String')))
set(handles.slider_tau, 'SliderStep', [SliderStep_k1/str2double(get(handles.edit_tau_max, 'String')) SliderStep_k2/str2double(get(handles.edit_tau_max, 'String'))])
if get(handles.slider_tau,'Value') < 0 
    set(handles.slider_tau,'Value', 0)
end
n_tau = get(handles.slider_tau,'Value');
set(handles.text5,'String', n_tau)

function SetAxes(handles)
%Применение настроек из переключателей (для осей)
if get(handles.radiobutton1, 'Value') == 1 %Заменить
	set(handles.axes1, 'NextPlot', 	'replace')
    set(handles.axes2, 'NextPlot', 	'replace')
    set(handles.axes3, 'NextPlot', 	'replace')
    set(handles.axes4, 'NextPlot', 	'replace')
    set(handles.axes1, 'NextPlot', 	'replace')
end
if get(handles.radiobutton1_1, 'Value') == 1 %Заменить в лимитах
	set(handles.axes1, 'NextPlot', 	'replacechildren')
    set(handles.axes2, 'NextPlot', 	'replacechildren')
    set(handles.axes3, 'NextPlot', 	'replacechildren')
    set(handles.axes4, 'NextPlot', 	'replacechildren')
end
if get(handles.radiobutton2, 'Value') == 1 %Добавить
	set(handles.axes1, 'NextPlot', 	'add')
    set(handles.axes2, 'NextPlot', 	'add')
    set(handles.axes3, 'NextPlot', 	'add')
    set(handles.axes4, 'NextPlot', 	'add')
end
if get(handles.radiobutton3, 'Value') == 1 %Лимиты осей
	axis(handles.axes1, 'auto'); %автоматические лимыты осей
    axis(handles.axes2, 'auto'); %автоматические лимыты осей
    axis(handles.axes3, 'auto'); %автоматические лимыты осей
    axis(handles.axes4, 'auto'); %автоматические лимыты осей
end
if get(handles.radiobutton4, 'Value') == 1 %Лимиты осей (жёстко заданные)
    achh_x_min = str2double(get(handles.edit_achh_x_min, 'String'));
    achh_x_max = str2double(get(handles.edit_achh_x_max, 'String'));
    achh_y_min = str2double(get(handles.edit_achh_y_min, 'String'));
    achh_y_max = str2double(get(handles.edit_achh_y_max, 'String'));
    axis(handles.axes1, [achh_x_min achh_x_max achh_y_min achh_y_max]);
    
    fchh_x_min = str2double(get(handles.edit_fchh_x_min, 'String'));
    fchh_x_max = str2double(get(handles.edit_fchh_x_max, 'String'));
    fchh_y_min = str2double(get(handles.edit_fchh_y_min, 'String'));
    fchh_y_max = str2double(get(handles.edit_fchh_y_max, 'String'));
    axis(handles.axes2, [fchh_x_min fchh_x_max fchh_y_min fchh_y_max]);
    
    afh_x_min = str2double(get(handles.edit_afh_x_min, 'String'));
    afh_x_max = str2double(get(handles.edit_afh_x_max, 'String'));
    afh_y_min = str2double(get(handles.edit_afh_y_min, 'String'));
    afh_y_max = str2double(get(handles.edit_afh_y_max, 'String'));
    axis(handles.axes3, [afh_x_min afh_x_max afh_y_min afh_y_max]);
    
    kr_x_min = str2double(get(handles.edit_kr_x_min, 'String'));
    kr_x_max = str2double(get(handles.edit_kr_x_max, 'String'));
    kr_y_min = str2double(get(handles.edit_kr_y_min, 'String'));
    kr_y_max = str2double(get(handles.edit_kr_y_max, 'String'));
    axis(handles.axes4, [kr_x_min kr_x_max kr_y_min kr_y_max]);
end

function SetOptions(handles)
SetSliders(handles)
SetAxes(handles)

function getresult(handles)
Draw(handles); %Отрисовка на графиках выбранного звена

% --- Executes during object creation, after setting all properties.
function slider_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%addlistener(hObject,'Value','PreSet',@(~,~)disp(get(hObject,'Value')));
%n = get(hObject,'Value');
%addlistener(hObject,'Value','PreSet',@(~,~)DoIt(n, handles));
%set(handles.text2,'String', n)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_k_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k_min as text
%        str2double(get(hObject,'String')) returns contents of edit_k_min as a double

% --- Executes during object creation, after setting all properties.
function edit_k_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k_max as text
%        str2double(get(hObject,'String')) returns contents of edit_k_max as a double

% --- Executes during object creation, after setting all properties.
function edit_k_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_z.
function popupmenu_z_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_z contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_z
global Num
%contents = get(handles.popupmenu_z,'String'); 
%Name = contents{get(handles.popupmenu_z, 'Value')}; % Получение названия
Num = get(handles.popupmenu_z, 'Value');
switch Num % свитч для видимости коэффициентов
    case 1
        set(handles.uipanel_k, 'Visible', 'on');
        set(handles.uipanel_t1_t2, 'Visible', 'off');
        set(handles.uipanel_ksi, 'Visible', 'off');
        set(handles.uipanel_t, 'Visible', 'off');
    case 2
        set(handles.uipanel_t1_t2, 'Visible', 'off');
        set(handles.uipanel_ksi, 'Visible', 'off');
        set(handles.uipanel_t, 'Visible', 'off');
        set(handles.uipanel_k, 'Visible', 'off');
    case 7
        set(handles.uipanel_t1_t2, 'Visible', 'off');
        set(handles.uipanel_ksi, 'Visible', 'off');
        set(handles.uipanel_t, 'Visible', 'off');
        set(handles.uipanel_k, 'Visible', 'on');
    case 10
        % Апериодич 2 порядка
        set(handles.uipanel_k, 'Visible', 'on');
        set(handles.uipanel_ksi, 'Visible', 'on');
        set(handles.uipanel_t1_t2, 'Visible', 'off');
        set(handles.uipanel_t, 'Visible', 'on');
    case 11
        % Интегро-дифференцирующее звено
        set(handles.uipanel_k, 'Visible', 'on');
        set(handles.uipanel_t1_t2, 'Visible', 'on');
        set(handles.uipanel_ksi, 'Visible', 'off');
        set(handles.uipanel_t, 'Visible', 'off');
    otherwise % Видимость по умолчанию
        set(handles.uipanel_k, 'Visible', 'on');
        set(handles.uipanel_t, 'Visible', 'on');
        set(handles.uipanel_ksi, 'Visible', 'off');
        set(handles.uipanel_t1_t2, 'Visible', 'off');
end
if  strcmp(get(handles.New_item_Clear_all, 'Checked'), 'on')
    Clear_All(handles)
end

% --- Executes during object creation, after setting all properties.
function popupmenu_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_k_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;

% --- Executes on slider movement.
function slider_tau_Callback(hObject, eventdata, handles)
% hObject    handle to slider_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles)
getresult(handles)

% --- Executes during object creation, after setting all properties.
function slider_tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_tau_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tau_min as text
%        str2double(get(hObject,'String')) returns contents of edit_tau_min as a double


% --- Executes during object creation, after setting all properties.
function edit_tau_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tau_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tau_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tau_max as text
%        str2double(get(hObject,'String')) returns contents of edit_tau_max as a double


% --- Executes during object creation, after setting all properties.
function edit_tau_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tau_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_ksi_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ksi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles)
getresult(handles)

% --- Executes during object creation, after setting all properties.
function slider_ksi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ksi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_ksi_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ksi_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ksi_min as text
%        str2double(get(hObject,'String')) returns contents of edit_ksi_min as a double


% --- Executes during object creation, after setting all properties.
function edit_ksi_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ksi_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ksi_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ksi_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ksi_max as text
%        str2double(get(hObject,'String')) returns contents of edit_ksi_max as a double


% --- Executes during object creation, after setting all properties.
function edit_ksi_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ksi_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_t_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles)
getresult(handles)

% --- Executes during object creation, after setting all properties.
function slider_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_t_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_min as text
%        str2double(get(hObject,'String')) returns contents of edit_t_min as a double


% --- Executes during object creation, after setting all properties.
function edit_t_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_t_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t_max as text
%        str2double(get(hObject,'String')) returns contents of edit_t_max as a double


% --- Executes during object creation, after setting all properties.
function edit_t_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4



function edit_achh_x_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_achh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_achh_x_min as text
%        str2double(get(hObject,'String')) returns contents of edit_achh_x_min as a double


% --- Executes during object creation, after setting all properties.
function edit_achh_x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_achh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_achh_x_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_achh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_achh_x_max as text
%        str2double(get(hObject,'String')) returns contents of edit_achh_x_max as a double


% --- Executes during object creation, after setting all properties.
function edit_achh_x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_achh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_achh_y_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_achh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_achh_y_max as text
%        str2double(get(hObject,'String')) returns contents of edit_achh_y_max as a double


% --- Executes during object creation, after setting all properties.
function edit_achh_y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_achh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_achh_y_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_achh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_achh_y_min as text
%        str2double(get(hObject,'String')) returns contents of edit_achh_y_min as a double


% --- Executes during object creation, after setting all properties.
function edit_achh_y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_achh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fchh_x_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fchh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fchh_x_min as text
%        str2double(get(hObject,'String')) returns contents of edit_fchh_x_min as a double


% --- Executes during object creation, after setting all properties.
function edit_fchh_x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fchh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fchh_x_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fchh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fchh_x_max as text
%        str2double(get(hObject,'String')) returns contents of edit_fchh_x_max as a double


% --- Executes during object creation, after setting all properties.
function edit_fchh_x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fchh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fchh_y_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fchh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fchh_y_max as text
%        str2double(get(hObject,'String')) returns contents of edit_fchh_y_max as a double


% --- Executes during object creation, after setting all properties.
function edit_fchh_y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fchh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fchh_y_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fchh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fchh_y_min as text
%        str2double(get(hObject,'String')) returns contents of edit_fchh_y_min as a double


% --- Executes during object creation, after setting all properties.
function edit_fchh_y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fchh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_afh_x_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_afh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_afh_x_min as text
%        str2double(get(hObject,'String')) returns contents of edit_afh_x_min as a double


% --- Executes during object creation, after setting all properties.
function edit_afh_x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_afh_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_afh_x_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_afh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_afh_x_max as text
%        str2double(get(hObject,'String')) returns contents of edit_afh_x_max as a double


% --- Executes during object creation, after setting all properties.
function edit_afh_x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_afh_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_afh_y_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_afh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_afh_y_max as text
%        str2double(get(hObject,'String')) returns contents of edit_afh_y_max as a double


% --- Executes during object creation, after setting all properties.
function edit_afh_y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_afh_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_afh_y_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_afh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_afh_y_min as text
%        str2double(get(hObject,'String')) returns contents of edit_afh_y_min as a double


% --- Executes during object creation, after setting all properties.
function edit_afh_y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_afh_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kr_x_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kr_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kr_x_min as text
%        str2double(get(hObject,'String')) returns contents of edit_kr_x_min as a double


% --- Executes during object creation, after setting all properties.
function edit_kr_x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kr_x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kr_x_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kr_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kr_x_max as text
%        str2double(get(hObject,'String')) returns contents of edit_kr_x_max as a double


% --- Executes during object creation, after setting all properties.
function edit_kr_x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kr_x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kr_y_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kr_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kr_y_max as text
%        str2double(get(hObject,'String')) returns contents of edit_kr_y_max as a double


% --- Executes during object creation, after setting all properties.
function edit_kr_y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kr_y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kr_y_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kr_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kr_y_min as text
%        str2double(get(hObject,'String')) returns contents of edit_kr_y_min as a double


% --- Executes during object creation, after setting all properties.
function edit_kr_y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kr_y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_slider_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slider_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slider_step as text
%        str2double(get(hObject,'String')) returns contents of edit_slider_step as a double


% --- Executes during object creation, after setting all properties.
function edit_slider_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slider_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_step_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step_time as text
%        str2double(get(hObject,'String')) returns contents of edit_step_time as a double


% --- Executes during object creation, after setting all properties.
function edit_step_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_kr_st_auto.
function checkbox_kr_st_auto_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_kr_st_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_kr_st_auto
SetOptions(handles);
getresult(handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function FindOut(handles)
if (get(handles.edit_achh_y_max,'String') == '2') && (get(handles.edit_achh_y_min,'String') == '0') && (get(handles.edit_fchh_y_min,'String') == '1') && (get(handles.edit_fchh_y_max,'String') == '6') 
    if (get(handles.radiobutton1_1,'Value') == 1)
        set(handles.About,'Visible', 'on');
    end
end

function edit_w_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_w_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_w_start as text
%        str2double(get(hObject,'String')) returns contents of edit_w_start as a double


% --- Executes during object creation, after setting all properties.
function edit_w_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_w_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_w_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_w_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_w_end as text
%        str2double(get(hObject,'String')) returns contents of edit_w_end as a double


% --- Executes during object creation, after setting all properties.
function edit_w_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_w_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_w_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_w_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_w_step as text
%        str2double(get(hObject,'String')) returns contents of edit_w_step as a double


% --- Executes during object creation, after setting all properties.
function edit_w_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_w_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function f_SaveToWorkspace
%assignin('base', 'A_base', a) 
%где a - переменная, созданная в функции, а A_base создается в рабочей среде.
prompt = {'Сохранить переменные с приставкой:'};
dlg_title = 'Приставка имён';
num_lines = 1;
defaultans = {''};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end

global  K Num T T1 T2 ksi tau w Wsis A Fi Ws Ws_p y t
assignin('base', [answer{1} 'K'], K);
if Num == 11 % если звено с T1 и T2
    assignin('base', [answer{1} 'T1'], T1);
    assignin('base', [answer{1} 'T2'], T2);
else % иначе все остальные звенья
    assignin('base', [answer{1} 'T'], T);
end
if Num == 10 % если звено с ksi
    assignin('base', [answer{1} 'ksi'], ksi);
end
assignin('base', [answer{1} 'tau'], tau);
assignin('base', [answer{1} 'w'], w);
assignin('base', [answer{1} 'Wsis'], Wsis);
assignin('base', [answer{1} 'A'], A);
assignin('base', [answer{1} 'Fi'], Fi);
assignin('base', [answer{1} 'Ws'], Ws);
assignin('base', [answer{1} 'Ws_p'], Ws_p);
assignin('base', [answer{1} 'y'], y);
assignin('base', [answer{1} 't'], t);

% --------------------------------------------------------------------
function SaveToWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to SaveToWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f_SaveToWorkspace;

% --------------------------------------------------------------------
function Grafik_Callback(hObject, eventdata, handles)
% hObject    handle to Grafik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ACHH_Callback(hObject, eventdata, handles)
% hObject    handle to ACHH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% A Fi Ws Ws_p
global Name Num w A K T T1 T2 ksi tau
title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
switch Num
    case 1 % если звено с К и tau
        title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 2 % если звено с tau
        title_text = [Name '. АЧХ. ' sprintf('tau=%1.2f', tau)];
    case 7 % если звено c К и tau
        title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 11 % если звено с T1 и T2
        title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, T1=%1.2f, T2=%1.2f tau=%1.2f', K, T1, T2, tau)];
    case 10 % если звено с ksi
        title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, T=%1.2f, ksi=%1.2f, tau=%1.2f', K, T, ksi, tau)];
    otherwise
        title_text = [Name '. АЧХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
end
figure('Name', title_text,'FileName',title_text);
plot(w, A, 'LineWidth', 2);
grid on
if  strcmp(get(handles.set_name, 'Checked'), 'on')
    title(title_text, 'FontSize', 15); 
end
if  strcmp(get(handles.set_axes, 'Checked'), 'on')
    xlabel('\omega', 'FontSize', 15); 
    ylabel('A(\omega)', 'Rotation', 0, 'FontSize', 15); 
end
if  strcmp(get(handles.set_legend, 'Checked'), 'on')
    switch Num % Подпись. Генерируется в MathType  с экспортом в буфер как Plain TeX
        case 1
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K$$');
            set(l,'Interpreter','latex')
        case 2
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {e^{ - p\tau }}$$');
            set(l,'Interpreter','latex')
        case 3
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 4
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 5
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {p \cdot \left( {T \cdot p + 1} \right)}}$$');
            set(l,'Interpreter','latex')
        case 6
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot \left( {T \cdot p + 1} \right)} \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 7
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot p$$');
            set(l,'Interpreter','latex')
        case 8
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot p} \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 9 
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot \left( {T \cdot p + 1} \right)$$');
            set(l,'Interpreter','latex')
        case 10
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {{T^2} \cdot {p^2} + 2 \cdot \zeta  \cdot T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 11
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot {{{T_1} \cdot p + 1} \over {{T_2} \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        otherwise
            l = legend('');
    end
    set(l, 'FontSize', 14);
end

% --------------------------------------------------------------------
function FCHH_Callback(hObject, eventdata, handles)
% hObject    handle to FCHH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% A Fi Ws Ws_p
global Name Num w Fi K T T1 T2 ksi tau
title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
switch Num
    case 1 % если звено с К и tau
        title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 2 % если звено с tau
        title_text = [Name '. ФЧХ. ' sprintf('tau=%1.2f', tau)];
    case 7 % если звено c К и tau
        title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 11 % если звено с T1 и T2
        title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, T1=%1.2f, T2=%1.2f tau=%1.2f', K, T1, T2, tau)];
    case 10 % если звено с ksi
        title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, T=%1.2f, ksi=%1.2f, tau=%1.2f', K, T, ksi, tau)];
    otherwise
        title_text = [Name '. ФЧХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
end
figure('Name', title_text,'FileName',title_text);
plot(w, Fi, 'LineWidth', 2);
grid on
if  strcmp(get(handles.set_name, 'Checked'), 'on')
    title(title_text, 'FontSize', 15); 
end
if  strcmp(get(handles.set_axes, 'Checked'), 'on')
    xlabel('\omega', 'FontSize', 15); 
    ylabel('\phi(\omega)', 'Rotation', 0, 'FontSize', 15); 
end
if  strcmp(get(handles.set_legend, 'Checked'), 'on')
    switch Num % Подпись. Генерируется в MathType  с экспортом в буфер как Plain TeX
        case 1
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K$$');
            set(l,'Interpreter','latex')
        case 2
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {e^{ - p\tau }}$$');
            set(l,'Interpreter','latex')
        case 3
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 4
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 5
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {p \cdot \left( {T \cdot p + 1} \right)}}$$');
            set(l,'Interpreter','latex')
        case 6
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot \left( {T \cdot p + 1} \right)} \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 7
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot p$$');
            set(l,'Interpreter','latex')
        case 8
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot p} \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 9 
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot \left( {T \cdot p + 1} \right)$$');
            set(l,'Interpreter','latex')
        case 10
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {{T^2} \cdot {p^2} + 2 \cdot \zeta  \cdot T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 11
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot {{{T_1} \cdot p + 1} \over {{T_2} \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        otherwise
            l = legend('');
    end
    set(l, 'FontSize', 14);
end

% --------------------------------------------------------------------
function AFH_Callback(hObject, eventdata, handles)
% hObject    handle to AFH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% A Fi Ws Ws_p
global Name Num K T T1 T2 ksi tau Ws Ws_p
title_text = [Name '. АФХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
switch Num
    case 1 % если звено с К и tau
        title_text = [Name '. АФХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 2 % если звено с tau
        title_text = [Name '. АФХ. ' sprintf('tau=%1.2f', tau)];
    case 7 % если звено c К и tau
        title_text = [Name '. АФХ. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 11 % если звено с T1 и T2
        title_text = [Name '. АФХ. ' sprintf('K=%1.2f, T1=%1.2f, T2=%1.2f, tau=%1.2f', K, T1, T2, tau)];
    case 10 % если звено с ksi
        title_text = [Name '. АФХ. ' sprintf('K=%1.2f, T=%1.2f, ksi=%1.2f, tau=%1.2f', K, T, ksi, tau)];
    otherwise
        title_text = [Name '. АФХ. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
end
figure('Name', title_text,'FileName',title_text);

plot(real(Ws), imag(Ws), 'LineWidth', 2);
if get(handles.checkbox_afh_ps, 'Value') == 1    %Показывать точки на АФХ
    hold 'on';
    plot(real(Ws), imag(Ws), '.', 'LineWidth', 2);
end
if get(handles.checkbox_afh_p, 'Value') == 1    %Начальная точка на АФХ
    hold 'on';
    plot(real(Ws_p), imag(Ws_p), 'Marker','o');
end
grid on
if  strcmp(get(handles.set_name, 'Checked'), 'on')
    title(title_text, 'FontSize', 15); 
end
if  strcmp(get(handles.set_axes, 'Checked'), 'on')
    xlabel('Re(\omega)', 'FontSize', 15); 
    ylabel('Im(\omega)', 'Rotation', 0, 'FontSize', 15);
end

if  strcmp(get(handles.set_legend, 'Checked'), 'on')
    switch Num % Подпись. Генерируется в MathType  с экспортом в буфер как Plain TeX
        case 1
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K$$');
            set(l,'Interpreter','latex')
        case 2
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {e^{ - p\tau }}$$');
            set(l,'Interpreter','latex')
        case 3
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 4
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 5
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {p \cdot \left( {T \cdot p + 1} \right)}}$$');
            set(l,'Interpreter','latex')
        case 6
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot \left( {T \cdot p + 1} \right)} \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 7
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot p$$');
            set(l,'Interpreter','latex')
        case 8
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot p} \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 9 
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot \left( {T \cdot p + 1} \right)$$');
            set(l,'Interpreter','latex')
        case 10
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {{T^2} \cdot {p^2} + 2 \cdot \zeta  \cdot T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 11
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot {{{T_1} \cdot p + 1} \over {{T_2} \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        otherwise
            l = legend('');
    end
    set(l, 'FontSize', 14);
end

% --------------------------------------------------------------------
function KR_Callback(hObject, eventdata, handles)
% hObject    handle to KR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% A Fi Ws Ws_p
global Name Num K T T1 T2 ksi tau y t
title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
switch Num
    case 1 % если звено с К и tau
        title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 2 % если звено с tau
        title_text = [Name '. Кривая разгона. ' sprintf('tau=%1.2f', tau)];
    case 7 % если звено c К и tau
        title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 11 % если звено с T1 и T2
        title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, T1=%1.2f, T2=%1.2f, tau=%1.2f', K, T1, T2, tau)];
    case 10 % если звено с ksi
        title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, T=%1.2f, ksi=%1.2f, tau=%1.2f', K, T, ksi, tau)];
    otherwise
        title_text = [Name '. Кривая разгона. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
end
figure('Name', title_text,'FileName', title_text);
plot( t, y, 'LineWidth', 2);
grid on
if  strcmp(get(handles.set_name, 'Checked'), 'on')
    title(title_text, 'FontSize', 15); 
end
if  strcmp(get(handles.set_axes, 'Checked'), 'on')
    xlabel('t', 'FontSize', 15); 
    ylabel('', 'Rotation', 0, 'FontSize', 15);  
end

if  strcmp(get(handles.set_legend, 'Checked'), 'on')
    switch Num % Подпись. Генерируется в MathType  с экспортом в буфер как Plain TeX
        case 1
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K$$');
            set(l,'Interpreter','latex')
        case 2
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {e^{ - p\tau }}$$');
            set(l,'Interpreter','latex')
        case 3
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 4
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 5
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {p \cdot \left( {T \cdot p + 1} \right)}}$$');
            set(l,'Interpreter','latex')
        case 6
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot \left( {T \cdot p + 1} \right)} \over {T \cdot p}}$$');
            set(l,'Interpreter','latex')
        case 7
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot p$$');
            set(l,'Interpreter','latex')
        case 8
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {{K \cdot p} \over {T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 9 
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot \left( {T \cdot p + 1} \right)$$');
            set(l,'Interpreter','latex')
        case 10
            l = legend('$${\mathop{\rm W}\nolimits} (p) = {K \over {{T^2} \cdot {p^2} + 2 \cdot \zeta  \cdot T \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        case 11
            l = legend('$${\mathop{\rm W}\nolimits} (p) = K \cdot {{{T_1} \cdot p + 1} \over {{T_2} \cdot p + 1}}$$');
            set(l,'Interpreter','latex')
        otherwise
            l = legend('');
    end
    set(l, 'FontSize', 14);
end

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function All_gr_Callback(hObject, eventdata, handles)
% hObject    handle to All_gr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Name Num K T T1 T2 ksi tau w A Fi Ws Ws_p y t
title_text = [Name '. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
switch Num
    case 1 % если звено с К и tau
        title_text = [Name '. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 2 % если звено с tau
        title_text = [Name '. ' sprintf('tau=%1.2f', tau)];
    case 7 % если звено c К и tau
        title_text = [Name '. ' sprintf('K=%1.2f, tau=%1.2f', K, tau)];
    case 11 % если звено с T1 и T2
        title_text = [Name '. ' sprintf('K=%1.2f, T1=%1.2f, T2=%1.2f tau=%1.2f', K, T1, T2, tau)];
    case 10 % если звено с ksi
        title_text = [Name '. ' sprintf('K=%1.2f, T=%1.2f, ksi=%1.2f, tau=%1.2f', K, T, ksi, tau)];
    otherwise
        title_text = [Name '. ' sprintf('K=%1.2f, T=%1.2f, tau=%1.2f', K, T, tau)];
end
figure('Name', title_text,'FileName',title_text);
subplot(2, 2, 1);
plot(w, A, 'LineWidth', 2);
grid on
title('АЧХ', 'FontSize', 15); 
xlabel('\omega', 'FontSize', 15); 
ylabel('A(\omega)', 'Rotation', 0, 'FontSize', 15); 

subplot(2, 2, 2);
plot(w, Fi, 'LineWidth', 2);
grid on
title('ФЧХ', 'FontSize', 15); 
xlabel('\omega', 'FontSize', 15); 
ylabel('\phi(\omega)', 'Rotation', 0, 'FontSize', 15); 

subplot(2, 2, 3);
plot(real(Ws), imag(Ws), 'LineWidth', 2);
hold 'on';
plot(real(Ws_p), imag(Ws_p), 'Marker','o');
grid on
title('АФХ', 'FontSize', 15); 
xlabel('Re(\omega)', 'FontSize', 15); 
ylabel('Im(\omega)', 'Rotation', 0, 'FontSize', 15);

subplot(2, 2, 4);
plot( t, y, 'LineWidth', 2);
grid on
title('Кривая разгона', 'FontSize', 15); 
xlabel('t', 'FontSize', 15); 
ylabel('', 'Rotation', 0, 'FontSize', 15);


% --------------------------------------------------------------------
function HelpFile_Callback(hObject, eventdata, handles)
% hObject    handle to HelpFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('Справка.docx');

% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'ТАУ от магистрантов - наглядное введение'; 'в Теорию Автоматического Управления!'; ''; 'Спасибо всем, кто принимал участие в разработке!'}, 'О программе');


% --- Executes on button press in checkbox_afh_p.
function checkbox_afh_p_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_afh_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_afh_p
SetOptions(handles);
getresult(handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetOptions(handles);
getresult(handles);


% --- Executes on key press with focus on pushbutton2 and none of its controls.
function pushbutton2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function edit_t1_t2_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t1_t2_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t1_t2_max as text
%        str2double(get(hObject,'String')) returns contents of edit_t1_t2_max as a double


% --- Executes during object creation, after setting all properties.
function edit_t1_t2_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t1_t2_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_t1_t2_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t1_t2_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t1_t2_min as text
%        str2double(get(hObject,'String')) returns contents of edit_t1_t2_min as a double


% --- Executes during object creation, after setting all properties.
function edit_t1_t2_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t1_t2_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_t1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles);
getresult(handles);

% --- Executes during object creation, after setting all properties.
function slider_t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_t2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    return
end
SetOptions(handles);
getresult(handles);

% --- Executes during object creation, after setting all properties.
function slider_t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function set_legend_Callback(hObject, eventdata, handles)
% hObject    handle to set_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'Checked'), 'on')
    set(hObject, 'Checked', 'off');
    set(handles.uitoggletool_legend, 'State', 'off');
else
    set(hObject, 'Checked', 'on');
    set(handles.uitoggletool_legend, 'State', 'on');
end


% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function New_item_Clear_all_Callback(hObject, eventdata, handles)
% hObject    handle to New_item_Clear_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(handles.New_item_Clear_all, 'Checked'), 'on')
    set(handles.New_item_Clear_all, 'Checked', 'off');
else
    set(handles.New_item_Clear_all, 'Checked', 'on');
end



function edit_step_k1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step_k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step_k1 as text
%        str2double(get(hObject,'String')) returns contents of edit_step_k1 as a double


% --- Executes during object creation, after setting all properties.
function edit_step_k1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step_k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_step_k2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step_k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step_k2 as text
%        str2double(get(hObject,'String')) returns contents of edit_step_k2 as a double


% --- Executes during object creation, after setting all properties.
function edit_step_k2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step_k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtool_save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f_SaveToWorkspace;


% --------------------------------------------------------------------
function uipushtool_get_result_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_get_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetOptions(handles);
getresult(handles);


% --------------------------------------------------------------------
function uipushtool_Clear_all_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_Clear_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Clear_All(handles);


% --------------------------------------------------------------------
function uitoggletool_legend_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'State'), 'on')
    set(hObject, 'State', 'on');
    set(handles.set_legend, 'Checked', 'on');
else
    set(hObject, 'State', 'off');
    set(handles.set_legend, 'Checked', 'off');
end


% --------------------------------------------------------------------
function set_name_Callback(hObject, eventdata, handles)
% hObject    handle to set_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'Checked'), 'on')
    set(hObject, 'Checked', 'off');
    set(handles.uitoggletool_name, 'State', 'off');
else
    set(hObject, 'Checked', 'on');
    set(handles.uitoggletool_name, 'State', 'on');
end

% --------------------------------------------------------------------
function set_axes_Callback(hObject, eventdata, handles)
% hObject    handle to set_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'Checked'), 'on')
    set(hObject, 'Checked', 'off');
    set(handles.uitoggletool_axes, 'State', 'off');
else
    set(hObject, 'Checked', 'on');
    set(handles.uitoggletool_axes, 'State', 'on');
end

% --------------------------------------------------------------------
function uitoggletool_name_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'State'), 'on')
    set(hObject, 'State', 'on');
    set(handles.set_name, 'Checked', 'on');
else
    set(hObject, 'State', 'off');
    set(handles.set_name, 'Checked', 'off');
end


% --------------------------------------------------------------------
function uitoggletool_axes_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(hObject, 'State'), 'on')
    set(hObject, 'State', 'on');
    set(handles.set_axes, 'Checked', 'on');
else
    set(hObject, 'State', 'off');
    set(handles.set_axes, 'Checked', 'off');
end


% --------------------------------------------------------------------
function uipushtool_achh_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_achh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ACHH_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool_FCHH_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_FCHH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FCHH_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool_AFH_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_AFH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AFH_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool_KR_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_KR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
KR_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uipushtool_All_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
All_gr_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox_afh_ps.
function checkbox_afh_ps_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_afh_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_afh_ps
SetOptions(handles);
getresult(handles);


% --------------------------------------------------------------------
function Control_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Control_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Zveno_Callback(hObject, eventdata, handles)
% hObject    handle to Zveno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function w_set_Callback(hObject, eventdata, handles)
% hObject    handle to w_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'С:','По:', 'Шаг:'};
dlg_title = 'Интервалы и шаг частоты';
num_lines = 1;
defaultans = {'0','20', '0.005'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
set(handles.edit_w_start, 'String', answer(1));
set(handles.edit_w_end, 'String', answer(2));
set(handles.edit_w_step, 'String', answer(3));

% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function k_set_Callback(hObject, eventdata, handles)
% hObject    handle to k_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'K:'};
dlg_title = 'K';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_k_min, 'String'));
    set(handles.edit_k_min, 'String', answer(1));
    set(handles.slider_k, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_k_max, 'String'));
    set(handles.edit_k_max, 'String', answer(1));
    set(handles.slider_k, 'Value', str2double(answer(1)));
else
    set(handles.slider_k, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function t_set_Callback(hObject, eventdata, handles)
% hObject    handle to t_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'T:'};
dlg_title = 'T';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_t_min, 'String'));
    set(handles.edit_t_min, 'String', answer(1));
    set(handles.slider_t, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_t_max, 'String'));
    set(handles.edit_t_max, 'String', answer(1));
    set(handles.slider_t, 'Value', str2double(answer(1)));
else
    set(handles.slider_t, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function tau_set_Callback(hObject, eventdata, handles)
% hObject    handle to tau_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'Тау:'};
dlg_title = 'Tау';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_tau_min, 'String'));
    set(handles.edit_tau_min, 'String', answer(1));
    set(handles.slider_tau, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_tau_max, 'String'));
    set(handles.edit_tau_max, 'String', answer(1));
    set(handles.slider_tau, 'Value', str2double(answer(1)));
else
    set(handles.slider_tau, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function ksi_set_Callback(hObject, eventdata, handles)
% hObject    handle to ksi_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'Кси:'};
dlg_title = 'Кси';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_ksi_min, 'String'));
    set(handles.edit_ksi_min, 'String', answer(1));
    set(handles.slider_ksi, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_ksi_max, 'String'));
    set(handles.edit_ksi_max, 'String', answer(1));
    set(handles.slider_ksi, 'Value', str2double(answer(1)));
else
    set(handles.slider_ksi, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function t1_set_Callback(hObject, eventdata, handles)
% hObject    handle to t1_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'T1:'};
dlg_title = 'T1';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_t1_t2_min, 'String'));
    set(handles.edit_t1_t2_min, 'String', answer(1));
    set(handles.slider_t1, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_t1_t2_max, 'String'));
    set(handles.edit_t1_t2_max, 'String', answer(1));
    set(handles.slider_t1, 'Value', str2double(answer(1)));
else
    set(handles.slider_t1, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function t2_set_Callback(hObject, eventdata, handles)
% hObject    handle to t2_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputdlg_w_set = {'T2:'};
dlg_title = 'T2';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(inputdlg_w_set,dlg_title,num_lines,defaultans);
if isempty(answer)
    return 
end
if str2double(answer(1)) < str2double(get(handles.edit_t1_t2_min, 'String'));
    set(handles.edit_t1_t2_min, 'String', answer(1));
    set(handles.slider_t2, 'Value', str2double(answer(1)));
else
if str2double(answer(1)) > str2double(get(handles.edit_t1_t2_max, 'String'));
    set(handles.edit_t1_t2_max, 'String', answer(1));
    set(handles.slider_t2, 'Value', str2double(answer(1)));
else
    set(handles.slider_t2, 'Value', str2double(answer(1)));
end
end
SetOptions(handles);
getresult(handles);

% --------------------------------------------------------------------
function w_set_start_Callback(hObject, eventdata, handles)
% hObject    handle to w_set_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Z1_Callback(hObject, eventdata, handles)
% hObject    handle to Z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 1);
popupmenu_z_Callback(hObject, eventdata, handles);

function slider_all_y_Move(handles)
uip_pos = get(handles.uipanel_all, 'Position');
dy = get(handles.slider_all_y,'Value');
set(handles.uipanel_all, 'Position', [uip_pos(1) 1.38-dy uip_pos(3) uip_pos(4)]);

function slider_all_x_Move(handles)
uip_pos = get(handles.uipanel_all, 'Position');
dx = get(handles.slider_all_x,'Value');
set(handles.uipanel_all, 'Position', [0-dx uip_pos(2) uip_pos(3) uip_pos(4)]);

% --- Executes on slider movement.
function slider_all_y_Callback(hObject, eventdata, handles)
% hObject    handle to slider_all_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_all_y_Move(handles)

% --- Executes during object creation, after setting all properties.
function slider_all_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_all_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function View_Callback(hObject, eventdata, handles)
% hObject    handle to View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_default_Callback(hObject, eventdata, handles)
% hObject    handle to view_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider_all_y, 'Value', 0.0);
set(handles.slider_all_x, 'Value', 0.0);
uip_pos = get(handles.uipanel_all, 'Position');
set(handles.uipanel_all, 'Position', [0 1.38 uip_pos(3) uip_pos(4)]);


% --------------------------------------------------------------------
function view_intervals_Callback(hObject, eventdata, handles)
% hObject    handle to view_intervals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_intervals_50_Callback(hObject, eventdata, handles)
% hObject    handle to view_intervals_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = 50;
set(handles.slider_all_y, 'Value', 0.0);
set(handles.slider_all_y, 'Min', -new_value);
set(handles.slider_all_y, 'Max', new_value);
set(handles.slider_all_y, 'SliderStep', [1/new_value 1])
set(handles.slider_all_x, 'Value', 0.0);
set(handles.slider_all_x, 'Min', -new_value);
set(handles.slider_all_x, 'Max', new_value);
set(handles.slider_all_x, 'SliderStep', [1/new_value 1])

% --------------------------------------------------------------------
function view_intervals_100_Callback(hObject, eventdata, handles)
% hObject    handle to view_intervals_100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = 100;
set(handles.slider_all_y, 'Value', 0.0);
set(handles.slider_all_y, 'Min', -new_value);
set(handles.slider_all_y, 'Max', new_value);
set(handles.slider_all_y, 'SliderStep', [1/new_value 1])
set(handles.slider_all_x, 'Value', 0.0);
set(handles.slider_all_x, 'Min', -new_value);
set(handles.slider_all_x, 'Max', new_value);
set(handles.slider_all_x, 'SliderStep', [1/new_value 1])

% --------------------------------------------------------------------
function view_intervals_200_Callback(hObject, eventdata, handles)
% hObject    handle to view_intervals_200 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = 200;
set(handles.slider_all_y, 'Value', 0.0);
set(handles.slider_all_y, 'Min', -new_value);
set(handles.slider_all_y, 'Max', new_value);
set(handles.slider_all_y, 'SliderStep', [1/new_value 1])
set(handles.slider_all_x, 'Value', 0.0);
set(handles.slider_all_x, 'Min', -new_value);
set(handles.slider_all_x, 'Max', new_value);
set(handles.slider_all_x, 'SliderStep', [1/new_value 1])

% --- Executes on slider movement.
function slider_all_x_Callback(hObject, eventdata, handles)
% hObject    handle to slider_all_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_all_x_Move(handles)

% --- Executes during object creation, after setting all properties.
function slider_all_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_all_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function Z2_Callback(hObject, eventdata, handles)
% hObject    handle to Z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 2);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z3_Callback(hObject, eventdata, handles)
% hObject    handle to Z3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 3);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z4_Callback(hObject, eventdata, handles)
% hObject    handle to Z4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 4);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z5_Callback(hObject, eventdata, handles)
% hObject    handle to Z5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 5);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z6_Callback(hObject, eventdata, handles)
% hObject    handle to Z6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 6);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z7_Callback(hObject, eventdata, handles)
% hObject    handle to Z7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 7);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z8_Callback(hObject, eventdata, handles)
% hObject    handle to Z8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 8);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z9_Callback(hObject, eventdata, handles)
% hObject    handle to Z9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 9);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z10_Callback(hObject, eventdata, handles)
% hObject    handle to Z10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 10);
popupmenu_z_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Z11_Callback(hObject, eventdata, handles)
% hObject    handle to Z11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu_z, 'Value', 11);
popupmenu_z_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function Disable_slider_result_Callback(hObject, eventdata, handles)
% hObject    handle to Disable_slider_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(handles.Disable_slider_result, 'Checked'), 'on')
    set(handles.Disable_slider_result, 'Checked', 'off');
else
    set(handles.Disable_slider_result, 'Checked', 'on');
end
