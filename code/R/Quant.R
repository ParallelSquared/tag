#   Q u a n t i t a t i o n   o f   D I A   -   a n a l y s i s   d u p l i c a t e d   f r o m   J M o d   m a n u s c r i p t 
 
 e v < -   r e a d . c s v ( " / V o l u m e s / L a b / K M D / J M o d _ T a g 6 P r e p r i n t / A l l _ R a t i o s . c s v " ) 
 
 h e a d ( e v ) 
 
 e v $ D a t a < - p a s t e 0 ( e v $ p l e x , " ,   " , e v $ A l g ) 
 
 g g p l o t ( e v [ e v $ C o m b i n a t i o n % i n % c ( " A / B " , " A / C " , " A / D " ) , ] , a e s ( y = R a t i o ,   x = C o m b i n a t i o n , f i l l = " X " ) )   + 
     g e o m _ b o x p l o t ( )   + 
     t h e m e _ t a g ( )   + 
     y l a b ( " L o g 2   Y e a s t   P r e c u r s o r   R a t i o s \ n " )   + 
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # f f 7 f 0 e " ) )   +   
     x l a b ( " " )   +   
     f a c e t _ w r a p ( " D a t a " , s h r i n k   =   T )   +   
     r r e m o v e ( " l e g e n d " )   +   
     y l i m ( c ( - 4 , 2 ) )   +   
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 4 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   + 
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 2 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   + 
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 3 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   #   +   
     # f o r c e _ p a n e l s i z e s ( c o l s   =   c ( 0 . 6 ,   1 ) )   
     
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / q u a n t _ d i a . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 1 . 3 ,   
               h e i g h t   =   i m a g e _ i n _ h * 1 . 3 , 
               u n i t s = " i n " ) 
 
 e v A D < - e v [ e v $ C o m b i n a t i o n % i n % c ( " A / D " ) , ] 
 
 g g p l o t ( e v A D , a e s ( y = R a t i o ,   x = C o m b i n a t i o n , f i l l = O r g ) )   + 
     g e o m _ b o x p l o t (   o u t l i e r . a l p h a   =   0 . 0 1 )   + 
     t h e m e _ t a g ( )   + 
     y l a b ( " L o g 2   Y e a s t   P r e c u r s o r   R a t i o s \ n " )   + 
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " g r e y 6 0 " , " # e 8 2 e 3 7 " ) )   +   
     x l a b ( " " )   +   
     f a c e t _ w r a p ( " D a t a " , s h r i n k   =   T )   +   
     r r e m o v e ( " l e g e n d . t i t l e " )   +   
     y l i m ( c ( - 4 , 4 ) )   +   
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 4 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   + 
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   #   +   
 # f o r c e _ p a n e l s i z e s ( c o l s   =   c ( 0 . 6 ,   1 ) )   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / q u a n t _ d i a 2 . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 
 g g p l o t ( e v [ e v $ C o m b i n a t i o n % i n % c ( " A / B " , " A / C " , " A / D " ) , ] , a e s ( y = R a t i o ,   x = C o m b i n a t i o n , f i l l = O r g ) )   + 
     g e o m _ b o x p l o t (   o u t l i e r . a l p h a   =   0 . 0 1 )   + 
     t h e m e _ t a g ( )   + 
     y l a b ( " L o g 2   Y e a s t   P r e c u r s o r   R a t i o s \ n " )   + 
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " g r e y 6 0 " , " # e 8 2 e 3 7 " ) )   +   
     x l a b ( " " )   +   
     f a c e t _ w r a p ( " D a t a " , s h r i n k   =   T )   +   
     r r e m o v e ( " l e g e n d . t i t l e " )   +   
     y l i m ( c ( - 4 , 4 ) )   +   
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 4 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   + 
     g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 3 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " ) +   #   +   
 g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 / 2 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " ) +   #   +   
 g e o m _ h l i n e ( y i n t e r c e p t   =   l o g 2 ( 1 ) ,   c o l = " g r e y 4 0 " ,   l t y = " d a s h e d " )   #   +   
 # f o r c e _ p a n e l s i z e s ( c o l s   =   c ( 0 . 6 ,   1 ) )   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / q u a n t _ d i a _ s u p p . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 2 . 6 ,   
               h e i g h t   =   i m a g e _ i n _ h * 1 . 3 , 
               u n i t s = " i n " ) 
 
 