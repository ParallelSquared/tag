#   S c r e e n i n g   c o m p o u n d s   o n   m o d e l   p e p t i d e s   f o r   i n c r e a s i n g   s p e c t r a l   q u a l i t y 
 
 #   S o u r c e   d e f a u l t   p a r a m e t e r s   f o r   a n a l y s i s 
 s o u r c e ( " p a r a m s . R " ) 
 
 #   L o a d   p a c k a g e s 
 l i b r a r y ( x m l 2 ) 
 l i b r a r y ( t i d y v e r s e ) 
 l i b r a r y ( g g p u b r ) 
 l i b r a r y ( g o o g l e d r i v e ) 
 l i b r a r y ( r e a d x l ) 
 
 #   R e a d   i n   t h e   d a t a b a s e   o f   m o d i f i c a t i o n s 
 c o m p o u n d s < - r e a d _ e x c e l ( " d a t a / M o d i f i c a t i o n s _ d b . x l s x " ) 
 
 r a w s < - c ( r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - A _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - B _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - C _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - D _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - E _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - F _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - G _ 1 0 n g " , 4 ) , 
                 r e p ( " 2 0 2 3 - 1 2 - 1 6 _ W a y n e _ 3 2 8 - 1 6 9 - H _ 1 0 n g " , 4 ) ) 
 
 c o r e s < - p a s t e 0 ( " T B B _ " , c ( " B " , " D " , " H " , " L " , 
                                               " E " , " I " , " M " , " P " , 
                                               " F " , " J " , " O " , " Q " , 
                                               " G " , " K " , " N " , " R " , 
                                               " B " , " D " , " H " , " L " , 
                                               " E " , " I " , " M " , " P " , 
                                               " F " , " J " , " O " , " Q " , 
                                               " G " , " K " , " N " , " R " ) ) 
 
 p e p s < - c ( r e p ( " A V Q V A H W A W S N E K " , 1 6 ) ,   r e p ( " T L E S S G E N L R L Q K " , 1 6 ) ) 
 
 c o p i e s < - c ( r e p ( 1 , 1 6 ) , r e p ( 2 , 1 6 ) ) 
 
 #   M e t a   d a t a   f o r   e x p e r i m e n t   1 6 9 
 m e t a 1 6 9 < - d a t a . f r a m e ( r a w s , c o r e s , p e p s , c o p i e s ) 
 
 #   P r e c i s i o n   o f   m a s s   m a t c h i n g 
 # d e c _ p t s < - 3 
 
 # e v 2 < - r e a d . d e l i m ( " / V o l u m e s / L a b / R e s u l t s / P e n t e l u t e / 2 0 2 4 - 0 1 - 1 9 _ a l l W a y n e / c o m b i n e d / t x t / m s m s S c a n s . t x t " ) 
 e v 2 < - r e a d . d e l i m ( " d a t a / e x p 1 6 9 _ t x t / m s m s S c a n s . t x t " ) 
 # e v 1 < - r e a d . d e l i m ( " / V o l u m e s / L a b / R e s u l t s / P e n t e l u t e / 2 0 2 4 - 0 1 - 1 9 _ a l l W a y n e / c o m b i n e d / t x t / a l l P e p t i d e s . t x t " ) 
 e v 1 < - r e a d . d e l i m ( " d a t a / e x p 1 6 9 _ t x t / a l l P e p t i d e s . t x t " ) 
 
 
 # # #   F i l t e r   t o   o n l y   E x p e r i m e n t   1 6 9 
 e v 2 < - e v 2 [ g r e p ( " 1 6 9 " , e v 2 $ R a w . f i l e ) ,   ] 
 e v 1 < - e v 1 [ g r e p ( " 1 6 9 " , e v 1 $ R a w . f i l e ) ,   ] 
 
 l f _ s c o r e _ p e p B < - m a x ( e v 1 $ S c o r e [ e v 1 $ S e q u e n c e % i n % " T L E S S G E N L R L Q K "   &   e v 1 $ D P . M o d i f i c a t i o n % i n % " U n m o d i f i e d " ] ) 
 l f _ s c o r e _ p e p A < - m a x ( e v 1 $ S c o r e [ e v 1 $ S e q u e n c e % i n % " A V Q V A H W A W S N E K "   &   e v 1 $ D P . M o d i f i c a t i o n % i n % " U n m o d i f i e d " ] ) 
 
 # #   c o m b i n e 
 
 e v 1 < - e v 1 [ ! i s . n a ( e v 1 $ D P . R a t i o . m o d . b a s e ) ,   ] 
 e v 1 < - e v 1 [ e v 1 $ D P . M a s s . D i f f e r e n c e > 1 0 0 ,   ] 
 
 
 e v 2 $ i d < - p a s t e 0 ( e v 2 $ R a w . f i l e , e v 2 $ D P . M o d . S c a n . N u m b e r ) 
 e v 1 $ i d < - p a s t e 0 ( e v 1 $ R a w . f i l e , e v 1 $ D P . M o d . S c a n . N u m b e r ) 
 
 e v < - m e r g e ( e v 1 , e v 2 , b y = " i d " ) 
 
 c o m p o u n d s < - c o m p o u n d s [ g r e p ( " T B B " , c o m p o u n d s $ N a m e ) ,   ] 
 c o m p o u n d s $ N u m P r i A m i d e < - c o m p o u n d s $ N u m P r i A m i d e - 1 
 c o m p o u n d s $ N u m S e c o n d A m i d e < - c o m p o u n d s $ N u m S e c o n d A m i d e + 1 
 
 e v [ , c o l n a m e s ( c o m p o u n d s ) ] < - N A 
 
 
 # # # # # #   m a t c h   r a w   f i l e s   a n d   p e p t i d e s 
 
 e v < - e v [ ( e v $ R a w . f i l e . x % i n % m e t a 1 6 9 $ r a w s [ 1 : 1 6 ]   &   e v $ D P . B a s e . S e q u e n c e . x % i n % m e t a 1 6 9 $ p e p s [ 1 : 1 6 ] ) | 
                   ( e v $ R a w . f i l e . x % i n % m e t a 1 6 9 $ r a w s [ 1 7 : 3 2 ]   &   e v $ D P . B a s e . S e q u e n c e . x % i n % m e t a 1 6 9 $ p e p s [ 1 7 : 3 2 ] ) , ] 
 
 e v $ n m o d < - 1 
 e v $ n m o d [ g r e p ( " T L E S S G E N L R L Q K " , e v $ D P . B a s e . S e q u e n c e . x ) ] < - 2 
 
 f o r ( i   i n   1 : n r o w ( c o m p o u n d s ) ) { 
     
     m i < - c o m p o u n d s $ M a s s [ i ] 
     
     m m < - e v $ D P . M a s s . D i f f e r e n c e . x   -   m i 
     
     i n d < - w h i c h ( a b s ( m m ) < ( m i / 1 e 5 ) )   #   m a t c h   a t   1 0 p p m   
     
     c i < - m e t a 1 6 9 $ c o r e s [ m e t a 1 6 9 $ r a w s % i n % e v $ R a w . f i l e . x [ i n d ] ] 
     
     m a t c h e s < - c ( ) 
     f o r ( j   i n   c i ) { 
         
         m a t c h e s < - c ( m a t c h e s ,   g r e p l ( j ,   c o m p o u n d s $ N a m e [ i ] ) ) 
         
     } 
     
     i f ( a n y ( m a t c h e s ) ) { 
         
         e v [ i n d ,   c o l n a m e s ( c o m p o u n d s ) ] < - c o m p o u n d s [ i , ] 
         
     } 
     
 } 
 
 e v f < - e v [ ! i s . n a ( e v $ M a s s ) ,   ] 
 
 c o l n a m e s ( e v f ) 
 
 e v f < - e v f [ , - 1 6 9 ] 
 
 
 
 e v f 2 < - e v f 
 
 
 e v a < - e v f [ e v f $ D P . B a s e . S e q u e n c e . y % i n % " A V Q V A H W A W S N E K " ,   ] 
 e v b < - e v f [ ! e v f $ D P . B a s e . S e q u e n c e . y % i n % " A V Q V A H W A W S N E K " ,   ] 
 
 
 p 1 < - g g p l o t ( e v b ,   a e s ( x = D P . S c o r e . y ) )   +   
     g e o m _ h i s t o g r a m ( b i n s   =   2 0 )   +   
     t h e m e _ t a g ( )   +   
     y l a b ( " C o u n t \ n " )   +   
     x l a b ( " \ n A n d r o m e d a   s c o r e " )   + 
     g e o m _ v l i n e ( x i n t e r c e p t   =   m a x ( e v b $ D P . S c o r e . y [ e v b $ N a m e % i n % " T a g   1 . 6 :   T B B _ J   +   B A 1 8 " ] ) ,   c o l o r = " r e d " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     # g e o m _ v l i n e ( x i n t e r c e p t   =   l f _ s c o r e _ p e p B ,   c o l o r = " g r e e n " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     r r e m o v e ( " l e g e n d . t i t l e " )   +   
     #   a n n o t a t e ( " t e x t " ,   x = 2 1 0 ,   y = 5 0 ,   l a b e l = " L e a d \ n c a n d i d a t e " ,   s i z e = 8 ,   c o l o r = " r e d " )   +   
     x l i m ( c ( 2 0 , 2 3 0 ) )   +   
     g g t i t l e ( " T L E S S G E N L R L Q K " ) 
 
 p 2 < - g g p l o t ( e v a ,   a e s ( x = D P . S c o r e . y ) )   +   
     g e o m _ h i s t o g r a m ( b i n s   =   2 0 )   +   
     t h e m e _ t a g ( )   +   
     y l a b ( " " )   +   
     x l a b ( " \ n A n d r o m e d a   s c o r e " )   + 
     g e o m _ v l i n e ( x i n t e r c e p t   =   m a x ( e v a $ D P . S c o r e . y [ e v a $ N a m e % i n % " T a g   1 . 6 :   T B B _ J   +   B A 1 8 " ] ) ,   c o l o r = " r e d " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     # g e o m _ v l i n e ( x i n t e r c e p t   =   l f _ s c o r e _ p e p A ,   c o l o r = " g r e e n " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     r r e m o v e ( " l e g e n d . t i t l e " )   +   
     # a n n o t a t e ( " t e x t " ,   x = 2 1 0 ,   y = 5 0 ,   l a b e l = " L e a d \ n c a n d i d a t e " ,   s i z e = 8 ,   c o l o r = " r e d " )   +   
     g g t i t l e ( " A V Q V A H W A W S N E K " ) 
 
 p 1   +   p 2   +   p a t c h w o r k : : p l o t _ l a y o u t ( n c o l = 2 ) 
 
 
 d f a < - e v a   % > %   g r o u p _ b y ( N a m e )   % > %   s u m m a r i s e ( m a x _ d p s c o r e   =   m a x ( D P . S c o r e . y ) ) 
 d f b < - e v b   % > %   g r o u p _ b y ( N a m e )   % > %   s u m m a r i s e ( m a x _ d p s c o r e   =   m a x ( D P . S c o r e . y ) ) 
 
 e v m < - m e r g e ( d f a , d f b , b y = " N a m e " ) 
 
 #   T a k e   u n i q u e ,   m a x 
 #   
 
 e v m $ c h o s e n < - " N o t   f u r t h e r   e v a l u a t e d " 
 e v m $ c h o s e n [ e v m $ N a m e % i n % " T a g   1 . 6 :   T B B _ J   +   B A 1 8 " ] < - " T o p   c h o i c e " 
 
 g g p l o t ( e v m ,   a e s ( x = m a x _ d p s c o r e . y ,   m a x _ d p s c o r e . x ,   c o l o r = c h o s e n ) )   +   
     g e o m _ p o i n t ( s i z e = 3 . 5 , a l p h a = 0 . 5 )   +   
     t h e m e _ t a g ( )   +   
     y l i m ( 0 , 3 0 0 )   +   
     x l i m ( 5 0 , 2 0 0 )   +   
     y l a b ( " " )   +   
     x l a b ( " \ n A n d r o m e d a   s c o r e \ n   T L E S S G E N L R L Q K " )   + 
     y l a b ( " A n d r o m e d a   s c o r e \ n   A V Q V A H W A W S N E K   \ n " )   + 
     s c a l e _ c o l o r _ m a n u a l ( v a l u e s = c ( " b l a c k " , " r e d " ) )   +   
     # g e o m _ v l i n e ( x i n t e r c e p t   =   m a x ( e v a $ D P . S c o r e . y [ e v a $ N a m e % i n % " T a g   1 . 6 :   T B B _ J   +   B A 1 8 " ] ) ,   c o l o r = " r e d " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     # g e o m _ v l i n e ( x i n t e r c e p t   =   l f _ s c o r e _ p e p A ,   c o l o r = " g r e e n " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     r r e m o v e ( " l e g e n d " )   
     # a n n o t a t e ( " t e x t " ,   x = 2 1 0 ,   y = 5 0 ,   l a b e l = " L e a d \ n c a n d i d a t e " ,   s i z e = 8 ,   c o l o r = " r e d " )   +   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / s c r e e n _ s c o r i n g 2 . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 
 e v l o g < - e v f [ e v f $ N a m e % i n % e v m $ N a m e , ] 
 e v l o g < - r e m o v e . d u p l i c a t e s ( e v l o g , " N a m e " ) 
 
 g g p l o t ( e v l o g ,   a e s ( x = L o g P ) )   +   
     g e o m _ h i s t o g r a m ( b i n s   =   2 0 )   +   
     t h e m e _ t a g ( )   +   
     y l a b ( " C o u n t \ n " )   +   
     x l a b ( " \ n L o g P " )   + 
     x l i m ( c ( - 2 , 9 ) )   + 
     y l i m ( c ( 0 , 2 0 ) )   +   
     t h e m e ( p a n e l . g r i d . m a j o r   =   e l e m e n t _ l i n e ( c o l o u r   =   " g r e y 9 0 " ,   s i z e   =   0 . 5 ) )   +   
     #   a n n o t a t e ( g e o m = " t e x t " ,   x = 5 ,   y = 1 5 ,   l a b e l = " T a g ,   L o g P   =   2 . 1 4 " , c o l o r = " r e d " , s i z e = 5 ,   h j u s t = 0 )   +   
     #   a n n o t a t e ( g e o m = " t e x t " ,   x = 5 ,   y = 1 3 ,   l a b e l = " G u a n i d i n e   d e r i v a t i v e ,   L o g P   =   0 . 9 3 " , c o l o r = " b l u e " , s i z e = 5 ,   h j u s t = 0 )   +   
     #   a n n o t a t e ( g e o m = " t e x t " ,   x = 5 ,   y = 1 1 ,   l a b e l = " T M T p r o ,   L o g P   =   0 . 0 6 " , c o l o r = " g r e y 2 0 " , s i z e = 5 ,   h j u s t = 0 )   +   
     #   a n n o t a t e ( g e o m = " t e x t " ,   x = 5 ,   y = 9 ,   l a b e l = " m T R A Q ,   L o g P   =   - 0 . 7 1 " , c o l o r = " g r e y 6 0 " , s i z e = 5 ,   h j u s t = 0 )   +   
     g e o m _ v l i n e ( x i n t e r c e p t   =   m e a n ( e v l o g [ e v l o g $ N a m e % i n % " T a g   1 . 6 :   T B B _ J   +   B A 1 8 " , " L o g P " ] ) ,   c o l o r = " r e d " ,   l t y = " d o t t e d " , s i z e = 1 . 5 )   +   
     g e o m _ v l i n e ( x i n t e r c e p t   =   - 0 . 7 1 ,   c o l o r = " g r e y 4 0 " ,   l t y = " d o t t e d " , s i z e = 1 . 5 )   +   
     g e o m _ v l i n e ( x i n t e r c e p t   =   1 . 4 2 ,   c o l o r = " b l u e " ,   l t y = " d o t t e d " , s i z e = 1 . 5 )   +   
     g e o m _ v l i n e ( x i n t e r c e p t   =   - 0 . 0 6 ,   c o l o r = " g r e y 2 0 " ,   l t y = " d o t t e d " , s i z e = 1 . 5 )   +   
     #   g e o m _ r e c t ( a e s ( x m i n   =   4 . 7 5 , 
     #                               x m a x   =   9 , 
     #                               y m i n   =   1 5 , 
     #                               y m a x   =   2 0 ) , 
     #                       f i l l   =   " w h i t e " ,   c o l o r   =   " g r e y 9 0 " ,   s i z e   =   1 )   +   
     a n n o t a t e ( g e o m = " t e x t " ,   x = 5 . 1 ,   y = 1 9 ,   l a b e l = " T a g " , c o l o r = " r e d " , s i z e = 4 ,   h j u s t = 0 )   +   
     a n n o t a t e ( g e o m = " t e x t " ,   x = 5 . 1 ,   y = 1 8 ,   l a b e l = " G u a n i d i n e   d e r i v a t i v e " , c o l o r = " b l u e " , s i z e = 4 ,   h j u s t = 0 )   +   
     a n n o t a t e ( g e o m = " t e x t " ,   x = 5 . 1 ,   y = 1 7 ,   l a b e l = " T M T p r o " , c o l o r = " g r e y 2 0 " , s i z e = 4 ,   h j u s t = 0 )   +   
     a n n o t a t e ( g e o m = " t e x t " ,   x = 5 . 1 ,   y = 1 6 ,   l a b e l = " m T R A Q " , c o l o r = " g r e y 4 0 " , s i z e = 4 ,   h j u s t = 0 )   +   
     # g e o m _ v l i n e ( x i n t e r c e p t   =   l f _ s c o r e _ p e p A ,   c o l o r = " g r e e n " ,   l t y = " d o t t e d " , s i z e = 2 )   +   
     r r e m o v e ( " l e g e n d . t i t l e " )   +   
     g e o m _ r e c t ( a e s ( x m i n   =   - 1 , 
                                 x m a x   =   1 , 
                                 y m i n   =   0 , 
                                 y m a x   =   2 0 ) , 
                         f i l l   =   " t r a n s p a r e n t " ,   c o l o r   =   " # 0 0 B A 3 8 " ,   s i z e   =   1 )   + 
     a n n o t a t e ( g e o m = " t e x t " ,   x = - 1 . 5 ,   y = 1 0 ,   l a b e l = " F u t u r e   c o m p o u n d   l i b r a r i e s " , c o l o r = " # 0 0 B A 3 8 " , s i z e = 4 ,   a n g l e = 9 0 )   
     
     # a n n o t a t e ( " t e x t " ,   x = 2 1 0 ,   y = 5 0 ,   l a b e l = " L e a d \ n c a n d i d a t e " ,   s i z e = 8 ,   c o l o r = " r e d " )   +   
     # g g t i t l e ( " A V Q V A H W A W S N E K " ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / l o g p 2 . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 