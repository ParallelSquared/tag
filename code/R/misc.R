# # #   M i s c e l l a n e o u s   a n a l y s i s   o f   T a g 
 
 #   L a d d e r   o f   D D M   l e v e l s   a n d   y i e l d   
 
 l i b r a r y ( a r r o w ) 
 l i b r a r y ( t i d y v e r s e ) 
 
 e v < - r e a d _ p a r q u e t ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / t a g / 2 0 0 p g / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 
 t a b l e ( e v $ R u n ) 
 
 e v $ D D M < - " 0 . 0 6 % " 
 e v $ D D M [ g r e p l ( " 1 2 D " , e v $ R u n ) ] < - " 0 . 1 2 % " 
 e v $ D D M [ g r e p l ( " 1 8 D " , e v $ R u n ) ] < - " 0 . 1 8 % " 
 e v $ D D M [ g r e p l ( " 2 4 D " , e v $ R u n ) ] < - " 0 . 2 4 % " 
 
 e v $ R e p < - " R e p   1 " 
 e v $ R e p [ g r e p l ( " r e p 2 " , e v $ R u n ) ] < - " R e p   2 " 
 
 g g p l o t ( e v , a e s ( x = D D M , f i l l = R e p ) )   +   
     g e o m _ b a r ( p o s i t i o n = " d o d g e " )   + 
     t h e m e _ t a g ( )   +   
     s c a l e _ f i l l _ g r e y ( )   +   
     r r e m o v e ( ' l e g e n d . t i t l e ' )   +   
     y l a b ( " P r e c u r s o r s ,   1 %   F D R \ n " )   +   
     x l a b ( " \ n   n - D o d e c y l - B e t a - M a l t o s i d e " ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / D D M _ l a d d e r _ I D s . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 # k p < - n a m e s ( t a b l e ( e v $ P r e c u r s o r . I d ) ) [ t a b l e ( e v $ P r e c u r s o r . I d ) = = 8 ] 
 # e v 2 < - e v [ e v $ P r e c u r s o r . I d % i n % k p ,   ]   % > %   g r o u p _ b y ( D D M ,   R e p )   % > %   s u m m a r i s e ( I n t   =   l o g 1 0 ( m e a n ( M s 1 . A r e a ) ) ) 
 
 e v 2 < - e v   % > %   p i v o t _ w i d e r ( i d _ c o l s   =   c ( P r e c u r s o r . I d ) ,   n a m e s _ f r o m   =   c ( D D M ) ,   v a l u e s _ f r o m   =   M s 1 . A r e a ,   v a l u e s _ f n   =   m e a n   ) 
 
 e v 2 [ , 2 ] < - e v 2 [ , 2 ] / e v 2 [ , 5 ] 
 e v 2 [ , 3 ] < - e v 2 [ , 3 ] / e v 2 [ , 5 ] 
 e v 2 [ , 4 ] < - e v 2 [ , 4 ] / e v 2 [ , 5 ] 
 
 e v 3 < - e v 2 [ , 1 : 4 ] 
 
 e v 3 < - e v 3   % > %   p i v o t _ l o n g e r ( ! P r e c u r s o r . I d ,   n a m e s _ t o   =   " D D M " ,   v a l u e s _ t o   =   " R a t i o " ) 
 
 e v 3 $ D D M [ e v 3 $ D D M % i n % " 0 . 1 2 % " ] < - " 0 . 1 2 % \ n   /   0 . 0 6 % " 
 e v 3 $ D D M [ e v 3 $ D D M % i n % " 0 . 1 8 % " ] < - " 0 . 1 8 % \ n   /   0 . 0 6 % " 
 e v 3 $ D D M [ e v 3 $ D D M % i n % " 0 . 2 4 % " ] < - " 0 . 2 4 % \ n   /   0 . 0 6 % " 
 
 g g p l o t ( e v 3 , a e s ( x = D D M ,   y = l o g 2 ( R a t i o ) ,   f i l l   =   D D M ) )   +   
     g e o m _ b o x p l o t ( o u t l i e r . s h a p e   =   N A ) + 
     y l i m ( c ( - 1 , 3 ) )   + 
     t h e m e _ t a g ( )   +   
   #   t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   3 0 ,   h j u s t = 0 . 5 ,   v j u s t = 0 . 5 ) )   +   
     s c a l e _ f i l l _ b r e w e r ( )   +   
     r r e m o v e ( ' l e g e n d ' )   +   
     y l a b ( e x p r e s s i o n ( ` R a t i o   o f   M S 1   A r e a ,   l o g ` [ 2 ] ) )   +   
     x l a b ( " \ n   n - D o d e c y l - B e t a - M a l t o s i d e " ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / D D M _ l a d d e r _ i n t . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 