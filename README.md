# Automatic-Control-Theory
Звенья теории автоматического управления (ТАУ).

Представлены звенья:  
Усилительное Звено 	tf(K, 'OutputDelay', tau)  
Звено чистого транспортного запаздывания 	tf(1, 'OutputDelay', tau)  
Апериодическое звено первого порядка 	tf(K, [T 1], 'OutputDelay', tau)  
Идеальное интегрирующее звено 	tf(K, [T 0], 'OutputDelay', tau)  
Реальное интегрирующее звено 	tf(K, [T 1 0], 'OutputDelay', tau)  
Пропорционально-интегральное звено 	tf([K*T K], [T 0], 'OutputDelay', tau)  
Идеальное дифференцирующее звено 	tf([K 0], 1, 'OutputDelay', tau)  
Реальное дифференцирующее звено (дифференцирующее инерционное звено) 	tf([K 0],[T 1],'OutputDelay', tau)  
Пропорционально-дифференциальное звено 	tf([K*T K], 1, 'OutputDelay', tau)  
Апериодическое звено второго порядка 	tf(K, [T^2 2*ksi*T 1], 'OutputDelay', tau)  
Интегро-дифференцирующее звено 	tf([K*T1 K], [T2 1], 'OutputDelay', tau)  

Video  
[![Automatic Control Theory, 1](https://i9.ytimg.com/vi/pbcMGc4qycg/mq2.jpg?sqp=CJiNvPkF&rs=AOn4CLAXBAj47x-vykUqlsAhlIzf3315bw)](https://youtu.be/pbcMGc4qycg "Automatic Control Theory, 1")  

[![Automatic Control Theory, 2](https://i9.ytimg.com/vi/qLkQUC7ibT0/mq2.jpg?sqp=CJiNvPkF&rs=AOn4CLA7LirsBJQKtFZi5095YWBdXIkNjA)](https://youtu.be/qLkQUC7ibT0 "Automatic Control Theory, 2")

Вся суть начинается с 120 строчки кода в файле TAU.m.  
Если выдаёт ошибку:  
Undefined function 'untitled' for input arguments of type 'char'.   
Error in @(hObject,eventdata)untitled('slider1_CreateFcn',hObject,eventdata,guidata(hObject))   
У меня тоже такая ошибка. Не знаю что это, такого элемента даже нет на форме.  
Но всё работает и с этой ошибкой.  

<a href="https://creativecommons.org/licenses/by/4.0/" Target="_blank"><img src="https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by.png" alt="CC BY" title="CC BY" width="120"></a>
