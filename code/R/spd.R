#   S a m p l e s   p e r   d a y 
 
 L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g - 5 p l e x " , " T a g - 5 p l e x " , " T a g - 9 p l e x " ) 
 G r a d i e n t < - c ( " 1 5   m i n . \ n g r a d i e n t " , " 3 0   m i n . \ n g r a d i e n t " , " 1 5   m i n . \ n g r a d i e n t " , " 3 0   m i n . \ n g r a d i e n t " , " 3 0   m i n . \ n g r a d i e n t " ) 
 S P D < - s i g n i f ( c ( 2 4 * 6 0 / 3 2 * 1 , 2 4 * 6 0 / 5 5 * 1 , 2 4 * 6 0 / 3 2 * 5 , 2 4 * 6 0 / 5 5 * 5 , 2 4 * 6 0 / 5 5 * 9   ) , 2 ) 
 T y p e < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , r e p ( " T a g " , 3 ) ) 
 d f < - d a t a . f r a m e ( L a b e l , G r a d i e n t , S P D , T y p e ) 
 
 d f $ L a b e l < - f a c t o r ( d f $ L a b e l ,   l e v e l s = c ( " L a b e l - f r e e " , " T a g - 5 p l e x " , " T a g - 9 p l e x " ) ) 
 
 g g p l o t ( d f , a e s ( x = L a b e l , y = S P D , f i l l = T y p e ) )   +   
     g e o m _ b a r ( s t a t   =   " i d e n t i t y " ,   p o s i t i o n   =   " d o d g e " , a e s ( f i l l = T y p e ) )   +   
     t h e m e _ t a g ( )   +   
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) )   +   
     x l a b ( " " )   +   
     y l a b ( " S a m p l e s   p e r   d a y   ( S P D ) \ n " )   +   
     g e o m _ t e x t ( a e s ( l a b e l = S P D , g r o u p = G r a d i e n t ) ,   p o s i t i o n = p o s i t i o n _ d o d g e ( w i d t h = 0 . 9 ) ,   v j u s t = - 0 . 5 , s i z e = 7 )   +   
     y l i m ( c ( 0 , 2 6 0 ) )   +   
     t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     t h e m e ( t e x t   =   e l e m e n t _ t e x t ( s i z e = 2 0 ) )   + 
     r r e m o v e ( " l e g e n d " )   + 
     f a c e t _ w r a p ( " G r a d i e n t " ,   s c a l e s   =   " f r e e _ x " )   +   
     f o r c e _ p a n e l s i z e s ( c o l s   =   c ( 2 ,   3 ) )   
     
     
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / s p d . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h * 1 . 5 , 
               u n i t s = " i n " ) 
 
 
 
 