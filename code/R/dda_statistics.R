#   D D A   s t a t i s t i c s 
 
 
 l f < - r e a d . d e l i m ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / a s t r a l _ C E _ D D A _ f i n a l / L F _ 2 8 / r e s u l t s . s a g e . t s v " ) 
 t 6 < - r e a d . d e l i m ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / a s t r a l _ C E _ D D A _ f i n a l / T 6 _ 2 4 / r e s u l t s . s a g e . t s v " ) 
 
 l f q < - r e a d . d e l i m ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / a s t r a l _ C E _ D D A _ f i n a l / L F _ 2 8 / l f q . t s v " ) 
 t 6 q < - r e a d . d e l i m ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / a s t r a l _ C E _ D D A _ f i n a l / T 6 _ 2 4 / l f q . t s v " ) 
 
 l f $ p e p t i d e < - p a s t e 0 ( l f $ p e p t i d e , l f $ c h a r g e ) 
 l f q $ p e p t i d e < - p a s t e 0 ( l f q $ p e p t i d e , l f q $ c h a r g e ) 
 t 6 $ p e p t i d e < - p a s t e 0 ( t 6 $ p e p t i d e , t 6 $ c h a r g e ) 
 t 6 q $ p e p t i d e < - p a s t e 0 ( t 6 q $ p e p t i d e , t 6 q $ c h a r g e ) 
 
 l f $ ` C - t e r m   a m i n o   a c i d ` < - " K " 
 l f $ ` C - t e r m   a m i n o   a c i d ` [ g r e p l ( " R " , l f $ p e p t i d e ) ] < - " R " 
 
 t 6 $ ` C - t e r m   a m i n o   a c i d ` < - " K " 
 t 6 $ ` C - t e r m   a m i n o   a c i d ` [ g r e p l ( " R " , t 6 $ p e p t i d e ) ] < - " R " 
 
 t 6 $ t y p e < - " T a g " 
 l f $ t y p e < - " L a b e l - f r e e " 
 
 f f 1 < - r b i n d ( t 6 , l f ) 
 
 
 d f _ t 6 < - d a t a . f r a m e ( t a b l e ( u n l i s t ( s t r _ s p l i t ( t 6 $ p e p t i d e , " " ) ) ) ) 
 d f _ l f < - d a t a . f r a m e ( t a b l e ( u n l i s t ( s t r _ s p l i t ( l f $ p e p t i d e , " " ) ) ) ) 
 
 d f < - m e r g e ( d f _ t 6 , d f _ l f ,   b y = " V a r 1 " ) 
 d f < - d f [ 1 4 : 3 3 , ] 
 d f $ F r e q . x < - d f $ F r e q . x / s u m ( d f $ F r e q . x ) 
 d f $ F r e q . y < - d f $ F r e q . y / s u m ( d f $ F r e q . y ) 
 c o l n a m e s ( d f ) < - c ( " V a r 1 " , " P S M t a g " , " L a b e l - f r e e " ) 
 
 d m < - m e l t ( d f ) 
 
 g g p l o t ( d m , a e s ( x = V a r 1 , y = v a l u e , f i l l = v a r i a b l e ) )   +   
     g e o m _ b a r ( s t a t = " i d e n t i t y " ,   p o s i t i o n = " d o d g e " )   +   
     t h e m e _ t a g ( )   +   
     g g t i t l e ( " R e l a t i v e   a m i n o   a c i d   f r e q u e n c i e s " )   +   
     x l a b ( " " )   +   
     y l a b ( " % " ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / d d a _ a a s . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 2 ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 g g p l o t ( f f 1 , a e s ( x = h y p e r s c o r e , a l p h a = ` C - t e r m   a m i n o   a c i d ` ,   f i l l = t y p e ) )   + 
     g e o m _ d e n s i t y ( )   +   
     s c a l e _ a l p h a _ d i s c r e t e ( r a n g e = c ( 0 . 7 , 0 . 2 ) )   +   
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) )   +   
     t h e m e _ t a g ( )   +   
     s c a l e _ x _ c o n t i n u o u s ( b r e a k s = s e q ( 1 0 , 5 5 , 5 ) )   +   
     x l i m ( c ( 1 5 , 5 5 ) ) +   
     # y l i m ( c ( 0 , 1 . 7 2 e 4 ) )   +   
     # t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     r r e m o v e ( " l e g e n d " )   +   
     # s c a l e _ x _ d i s c r e t e ( l a b e l s = c ( " P S M s " , " I d e n t i f i e d   p e p t i d e s " , " Q u a n t i f i e d   p e p t i d e s " , " I d e n t i f i e d   p r o t e i n s " , " Q u a n t i f i e d   p r o t e i n s " ) )   +   
     x l a b ( " \ n H y p e r s c o r e " )   + 
     y l a b ( " D e n s i t y \ n " )   +   
     # s c a l e _ f i l l _ g r e y ( )   +   
     f a c e t _ w r a p ( " t y p e " ,   s c a l e s   =   " f r e e " )   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / d d a _ s c o r e s . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 2 ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 l f f < - l f [ l f $ s p e c t r u m _ q < 0 . 0 1 , ] 
 t 6 f < - t 6 [ t 6 $ s p e c t r u m _ q < 0 . 0 1 , ] 
 
 #   l f f < - l f [ l f $ p r o t e i n _ q < 0 . 0 1 , ] 
 #   t 6 f < - t 6 [ t 6 $ p r o t e i n _ q < 0 . 0 1 , ] 
 
 
 l f f 1 < - l f f 
 t 6 f 1 < - t 6 f 
 
 #   K 
 l f f < - l f f 1 [ g r e p l ( " K " , l f f 1 $ p e p t i d e ) , ] 
 t 6 f < - t 6 f 1 [ g r e p l ( " K " , t 6 f 1 $ p e p t i d e ) , ] 
 
 i d e n t i f i e d _ p e p t i d e s < - c ( l e n g t h ( u n i q u e ( l f f $ p e p t i d e ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p e p t i d e ) ) ) 
 i d e n t i f i e d _ p r o t < - c ( l e n g t h ( u n i q u e ( l f f $ p r o t e i n s ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p r o t e i n s ) ) ) 
 p s m s < - c ( n r o w ( l f f ) , n r o w ( t 6 f ) ) 
 q u a n t i f i e d _ p e p t i d e s < - c ( l e n g t h ( u n i q u e ( l f f $ p e p t i d e [ l f f $ p e p t i d e % i n % l f q $ p e p t i d e ] ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p e p t i d e [ t 6 f $ p e p t i d e % i n % t 6 q $ p e p t i d e ] ) ) ) 
 q u a n t i f i e d _ p r o t < - c ( l e n g t h ( u n i q u e ( l f f $ p r o t e i n s [ l f f $ p e p t i d e % i n % l f q $ p e p t i d e ] ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p r o t e i n s [ t 6 f $ p e p t i d e % i n % t 6 q $ p e p t i d e ] ) ) ) 
 t y p e < - c ( " L a b e l - f r e e " , " T a g " ) 
 a a < - c ( " K " , " K " ) 
 
 # d f < - d a t a . f r a m e ( p s m s ,   i d e n t i f i e d _ p e p t i d e s ,   i d e n t i f i e d _ p r o t , q u a n t i f i e d _ p e p t i d e s ,   q u a n t i f i e d _ p r o t , t y p e , a a ) ;   d f 
 d f < - d a t a . f r a m e ( p s m s ,   i d e n t i f i e d _ p e p t i d e s ,   i d e n t i f i e d _ p r o t , t y p e , a a ) ;   d f 
 
 
 #   R 
 l f f < - l f f 1 [ g r e p l ( " R " , l f f 1 $ p e p t i d e ) , ] 
 t 6 f < - t 6 f 1 [ g r e p l ( " R " , t 6 f 1 $ p e p t i d e ) , ] 
 
 i d e n t i f i e d _ p e p t i d e s < - c ( l e n g t h ( u n i q u e ( l f f $ p e p t i d e ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p e p t i d e ) ) ) 
 i d e n t i f i e d _ p r o t < - c ( l e n g t h ( u n i q u e ( l f f $ p r o t e i n s ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p r o t e i n s ) ) ) 
 p s m s < - c ( n r o w ( l f f ) , n r o w ( t 6 f ) ) 
 q u a n t i f i e d _ p e p t i d e s < - c ( l e n g t h ( u n i q u e ( l f f $ p e p t i d e [ l f f $ p e p t i d e % i n % l f q $ p e p t i d e ] ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p e p t i d e [ t 6 f $ p e p t i d e % i n % t 6 q $ p e p t i d e ] ) ) ) 
 q u a n t i f i e d _ p r o t < - c ( l e n g t h ( u n i q u e ( l f f $ p r o t e i n s [ l f f $ p e p t i d e % i n % l f q $ p e p t i d e ] ) ) ,   l e n g t h ( u n i q u e ( t 6 f $ p r o t e i n s [ t 6 f $ p e p t i d e % i n % t 6 q $ p e p t i d e ] ) ) ) 
 t y p e < - c ( " L a b e l - f r e e " , " T a g " ) 
 a a < - c ( " R " , " R " ) 
 
 # d f 2 < - d a t a . f r a m e ( p s m s ,   i d e n t i f i e d _ p e p t i d e s ,   i d e n t i f i e d _ p r o t , q u a n t i f i e d _ p e p t i d e s ,   q u a n t i f i e d _ p r o t , t y p e , a a ) ;   d f 
 d f 2 < - d a t a . f r a m e ( p s m s ,   i d e n t i f i e d _ p e p t i d e s ,   i d e n t i f i e d _ p r o t , t y p e , a a ) ;   d f 
 
 d f f < - r b i n d ( d f , d f 2 ) 
 
 l i b r a r y ( r e s h a p e 2 ) 
 d f m < - m e l t ( d f f ) 
 
 c o l n a m e s ( d f m ) [ 2 ] < - " C - t e r m   a m i n o   a c i d " 
 
 # d f m $ v a r i a b l e < - f a c t o r ( d f m $ v a r i a b l e ,   l e v e l s = c ( " p s m s " , " i d e n t i f i e d _ p e p t i d e s " , " q u a n t i f i e d _ p e p t i d e s " , " i d e n t i f i e d _ p r o t " , " q u a n t i f i e d _ p r o t " ) ) 
 d f m $ v a r i a b l e < - f a c t o r ( d f m $ v a r i a b l e ,   l e v e l s = c ( " p s m s " , " i d e n t i f i e d _ p e p t i d e s " , " i d e n t i f i e d _ p r o t " ) ) 
 
 g g p l o t ( d f m , a e s ( x = v a r i a b l e , y = v a l u e , a l p h a = ` C - t e r m   a m i n o   a c i d ` ,   f i l l = t y p e ) )   + 
     g e o m _ b a r ( s t a t = " i d e n t i t y " ,   p o s i t i o n = " d o d g e " )   + 
     s c a l e _ a l p h a _ d i s c r e t e ( r a n g e = c ( 0 . 7 , 0 . 3 ) )   +   
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) )   +   
     t h e m e _ t a g ( )   +   
     # y l i m ( c ( 0 , 1 . 7 2 e 4 ) )   +   
     t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     # r r e m o v e ( " l e g e n d " )   +   
     g u i d e s ( f i l l = " n o n e " )   + 
   #   s c a l e _ x _ d i s c r e t e ( l a b e l s = c ( " P S M s " , " I d e n t i f i e d   p e p t i d e s " , " Q u a n t i f i e d   p e p t i d e s " , " I d e n t i f i e d   p r o t e i n s " , " Q u a n t i f i e d   p r o t e i n s " ) )   +   
     s c a l e _ x _ d i s c r e t e ( l a b e l s = c ( " P S M s " , " P e p t i d e s " , " P r o t e i n s " ) )   +   
     x l a b ( " " )   + 
     y l a b ( " C o u n t ,   1   % F D R \ n " )   +   
     # s c a l e _ f i l l _ g r e y ( )   +   
     f a c e t _ w r a p ( ~ t y p e ,   n c o l = 2 )   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / d d a _ i d s . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 1 ,   
               h e i g h t   =   i m a g e _ i n _ h * 2 , 
               u n i t s = " i n " ) 
 
 t 6 f 1 $ t y p e < - " T a g " 
 l f f 1 $ t y p e < - " L a b e l - f r e e " 
 
 f f < - r b i n d ( l f f 1 , t 6 f 1 ) 
 
 g g p l o t ( f f , a e s ( x = r t , a l p h a = ` C - t e r m   a m i n o   a c i d ` ,   f i l l = t y p e ) )   + 
     g e o m _ d e n s i t y ( )   +   
     s c a l e _ a l p h a _ d i s c r e t e ( r a n g e = c ( 0 . 7 , 0 . 2 ) )   +   
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) )   +   
     t h e m e _ t a g ( )   +   
     s c a l e _ x _ c o n t i n u o u s ( b r e a k s = s e q ( 1 0 , 5 5 , 5 ) )   +   
     x l i m ( c ( 1 5 , 5 5 ) ) +   
     # y l i m ( c ( 0 , 1 . 7 2 e 4 ) )   +   
     # t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     r r e m o v e ( " l e g e n d " )   +   
     # s c a l e _ x _ d i s c r e t e ( l a b e l s = c ( " P S M s " , " I d e n t i f i e d   p e p t i d e s " , " Q u a n t i f i e d   p e p t i d e s " , " I d e n t i f i e d   p r o t e i n s " , " Q u a n t i f i e d   p r o t e i n s " ) )   +   
     x l a b ( " \ n R e t e n t i o n   t i m e   ( m i n u t e s ) " )   + 
     y l a b ( " D e n s i t y ,   1 %   F D R \ n " )   +   
     # s c a l e _ f i l l _ g r e y ( )   +   
     f a c e t _ w r a p ( " t y p e " ,   s c a l e s   =   " f r e e " )   
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / d d a _ r t . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w * 2 ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 #   f f $ p e p t i d e < - g s u b ( " \ \ [ ( . * ) \ \ ] + " , " " ,   f f $ p e p t i d e ) 
 #   f f $ p e p t i d e < - g s u b ( " ( [ 0 - 9 ] ) " , " " ,   f f $ p e p t i d e ) 
 
 f f $ p e p t i d e < - g s u b ( " [ + 3 0 8 . 1 1 6 1 ] - " , " " , f f $ p e p t i d e , f i x e d = T ) 
 f f $ p e p t i d e < - g s u b ( " [ + 3 0 8 . 1 1 6 1 ] " , " " , f f $ p e p t i d e , f i x e d = T ) 
 f f $ p e p t i d e < - g s u b ( " ( [ 0 - 9 ] ) " , " " ,   f f $ p e p t i d e ) 
 
 
 #   t o t p e p < - u n i q u e ( f f $ p e p t i d e ) 
 #   t o t p r o t < - u n i q u e ( f f $ p r o t e i n s ) 
 
 # i n s t a l l . p a c k a g e s ( " g g V e n n D i a g r a m " ) 
 l i b r a r y ( g g V e n n D i a g r a m ) 
 
 
 l f 2 4 < - r e a d . d e l i m ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / p r o c e s s e d _ d a t a / a s t r a l _ C E _ D D A _ f i n a l / L F _ 2 4 / r e s u l t s . s a g e . t s v " ) 
 l f 2 4 < - l f 2 4 [ l f 2 4 $ s p e c t r u m _ q < 0 . 0 1 , ] 
 
 #   x < - l i s t ( 
 #       L F = f f $ p e p t i d e [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] , 
 #       L F = l f 2 4 $ p e p t i d e 
 #   ) 
 #   
 #   g g V e n n D i a g r a m ( x ,   s e t _ s i z e   =   8 ,   l a b e l _ s i z e   =   8 )   +   g g t i t l e ( " P e p t i d e s ,   1 %   F D R " )   + t h e m e ( t e x t = e l e m e n t _ t e x t ( s i z e = 2 5 ) )   + 
 #       s c a l e _ f i l l _ g r a d i e n t n ( c o l o u r s = g r e y . c o l o r s ( 2 ,   s t a r t   =   1 ,   e n d   = 0 . 5 ) , l i m i t s = c ( 0 , 1 0 0 0 0 ) )   +   
 #       r r e m o v e ( " l e g e n d " )   +       s c a l e _ x _ c o n t i n u o u s ( e x p a n d   =   e x p a n s i o n ( m u l t   =   . 1 ) ) 
 #   
 #   #   g g V e n n D i a g r a m ( x ,   s e t _ s i z e   =   8 ,   l a b e l _ s i z e   =   8 )   +   g g t i t l e ( " P e p t i d e s ,   1 %   F D R " )   + t h e m e ( t e x t = e l e m e n t _ t e x t ( s i z e = 2 5 ) )   + 
 #   #       # s c a l e _ f i l l _ g r a d i e n t n ( c o l o u r s = g r e y . c o l o r s ( 2 ,   s t a r t   =   1 ,   e n d   = 0 . 5 ) , l i m i t s = c ( 0 , 1 0 0 0 0 ) )   +   
 #   #       s c a l e _ f i l l _ g r a d i e n t 2 ( l o w   =   , l i m i t s = c ( 0 , 1 0 0 0 0 ) )   +   
 #   #       r r e m o v e ( " l e g e n d " )   +       s c a l e _ x _ c o n t i n u o u s ( e x p a n d   =   e x p a n s i o n ( m u l t   =   . 1 ) ) 
 #   
 #   
 #   g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n 1 . p d f " ,   
 #                 d e v i c e = " p d f " ,   
 #                 w i d t h   =   i m a g e _ i n _ w ,   
 #                 h e i g h t   =   i m a g e _ i n _ h , 
 #                 u n i t s = " i n " ) 
 #   
 #   x < - l i s t ( 
 #       L F = f f $ p r o t e i n s [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] , 
 #       L F = l f 2 4 $ p r o t e i n s 
 #   ) 
 #   
 #   g g V e n n D i a g r a m ( x ,   s e t _ s i z e   =   8 ,   l a b e l _ s i z e   =   8 )   +   g g t i t l e ( " P r o t e i n s ,   1 %   F D R " )   + t h e m e ( t e x t = e l e m e n t _ t e x t ( s i z e = 2 5 ) )   + 
 #       s c a l e _ f i l l _ g r a d i e n t n ( c o l o u r s = g r e y . c o l o r s ( 2 ,   s t a r t   =   1 ,   e n d   = 0 . 5 ) , l i m i t s = c ( 0 , 4 0 0 0 ) ) +   r r e m o v e ( " l e g e n d " )   +       s c a l e _ x _ c o n t i n u o u s ( e x p a n d   =   e x p a n s i o n ( m u l t   =   . 1 ) ) 
 #   
 #   g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n 2 . p d f " ,   
 #                 d e v i c e = " p d f " ,   
 #                 w i d t h   =   i m a g e _ i n _ w ,   
 #                 h e i g h t   =   i m a g e _ i n _ h , 
 #                 u n i t s = " i n " ) 
 #   
 #   
 #   x < - l i s t ( 
 #       L F = f f $ p e p t i d e [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] , 
 #       T a g = f f $ p e p t i d e [ ! f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
 #   ) 
 #   
 #   g g V e n n D i a g r a m ( x ,   s e t _ s i z e   = 8   ,   l a b e l _ s i z e   =   8 )   +   g g t i t l e ( " P e p t i d e s ,   1 %   F D R " )   + t h e m e ( t e x t = e l e m e n t _ t e x t ( s i z e = 2 5 ) )   + 
 #   s c a l e _ f i l l _ g r a d i e n t n ( c o l o u r s = g r e y . c o l o r s ( 2 ,   s t a r t   =   1 ,   e n d   = 0 . 5 ) , l i m i t s = c ( 0 , 1 0 0 0 0 ) ) +   r r e m o v e ( " l e g e n d " )   +   
 #       s c a l e _ x _ c o n t i n u o u s ( e x p a n d   =   e x p a n s i o n ( m u l t   =   . 1 ) ) 
 #   
 #   g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n 3 . p d f " ,   
 #                 d e v i c e = " p d f " ,   
 #                 w i d t h   =   i m a g e _ i n _ w ,   
 #                 h e i g h t   =   i m a g e _ i n _ h , 
 #                 u n i t s = " i n " ) 
 #   
 #   x < - l i s t ( 
 #       L F = f f $ p r o t e i n s [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] , 
 #       T a g = f f $ p r o t e i n s [ ! f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
 #   ) 
 #   
 #   g g V e n n D i a g r a m ( x ,   s e t _ s i z e   =   8 ,   l a b e l _ s i z e   =   8 )   +   g g t i t l e ( " P r o t e i n s ,   1 %   F D R " )   + t h e m e ( t e x t = e l e m e n t _ t e x t ( s i z e = 2 5 ) )   +   
 #       s c a l e _ f i l l _ g r a d i e n t n ( c o l o u r s = g r e y . c o l o r s ( 2 ,   s t a r t   =   1 ,   e n d   = 0 . 5 ) , l i m i t s = c ( 0 , 4 0 0 0 ) ) +   r r e m o v e ( " l e g e n d " )   + 
 #       s c a l e _ x _ c o n t i n u o u s ( e x p a n d   =   e x p a n s i o n ( m u l t   =   . 1 ) ) 
 #   
 #   
 #   g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n 4 . p d f " ,   
 #                 d e v i c e = " p d f " ,   
 #                 w i d t h   =   i m a g e _ i n _ w ,   
 #                 h e i g h t   =   i m a g e _ i n _ h , 
 #                 u n i t s = " i n " ) 
 
 
 
 # # # #   V e n n   1 
     L F 1 = f f $ p e p t i d e [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
     L F 2 = l f 2 4 $ p e p t i d e 
     
 T L = u n i q u e ( i n t e r s e c t ( L F 1 , L F 2 ) ) 
 l f 1 < - u n i q u e ( s e t d i f f ( L F 1 , L F 2 ) ) 
 l f 2 < - u n i q u e ( s e t d i f f ( L F 2 , L F 1 ) ) 
 
 p 1 < - p l o t ( e u l e r ( c ( " L F "   =   l e n g t h ( l f 1 ) ,   " L F   "   =   l e n g t h ( l f 2 ) ,   " L F & L F   "   =   l e n g t h ( T L ) ) ,   s h a p e   =   " e l l i p s e " ) ,   
                   f i l l = c ( " # 1 f 7 7 b 4 " , " # 1 f 7 7 b 4 " ) , 
                   a l p h a = 0 . 6 , 
                   q u a n t i t i e s   =   l i s t ( c e x = 2 ) , 
                   l a b e l s   =   l i s t ( c e x = 2 ) , 
                   t i t l e = " " ) 
 
 g < - g r i d . a r r a n g e ( g r o b s   =   l i s t ( p 1 ) ,   t o p   =   t e x t G r o b ( " P e p t i d e s ,   1 %   F D R " ,   g p = g p a r ( f o n t s i z e = 2 5 , f o n t = 8 ) ) ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n - 1 . p d f " , 
               g , 
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 # # # # 
 
 # # # #   V e n n   2 
     L F 1 = f f $ p r o t e i n s [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
     L F 2 = l f 2 4 $ p r o t e i n s 
 
 T L = u n i q u e ( i n t e r s e c t ( L F 1 , L F 2 ) ) 
 l f 1 < - u n i q u e ( s e t d i f f ( L F 1 , L F 2 ) ) 
 l f 2 < - u n i q u e ( s e t d i f f ( L F 2 , L F 1 ) ) 
 
 p 1 < - p l o t ( e u l e r ( c ( " L F "   =   l e n g t h ( l f 1 ) ,   " L F   "   =   l e n g t h ( l f 2 ) ,   " L F & L F   "   =   l e n g t h ( T L ) ) ,   s h a p e   =   " e l l i p s e " ) ,   
                   f i l l = c ( " # 1 f 7 7 b 4 " , " # 1 f 7 7 b 4 " ) , 
                   a l p h a = 0 . 6 , 
                   q u a n t i t i e s   =   l i s t ( c e x = 2 ) , 
                   l a b e l s   =   l i s t ( c e x = 2 ) , 
                   t i t l e = " " ) 
 
 g < - g r i d . a r r a n g e ( g r o b s   =   l i s t ( p 1 ) ,   t o p   =   t e x t G r o b ( " P r o t e i n s ,   1 %   F D R " ,   g p = g p a r ( f o n t s i z e = 2 5 , f o n t = 8 ) ) ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n - 2 . p d f " , 
               g , 
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 # # # # 
 
 
 
 # # # #   V e n n   3 
     L F 1 = f f $ p e p t i d e [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
     L F 2 = f f $ p e p t i d e [ ! f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
 
 T L = u n i q u e ( i n t e r s e c t ( L F 1 , L F 2 ) ) 
 l f 1 < - u n i q u e ( s e t d i f f ( L F 1 , L F 2 ) ) 
 l f 2 < - u n i q u e ( s e t d i f f ( L F 2 , L F 1 ) ) 
 
 p 1 < - p l o t ( e u l e r ( c ( " L F "   =   l e n g t h ( l f 1 ) ,   " T a g "   =   l e n g t h ( l f 2 ) ,   " L F & T a g "   =   l e n g t h ( T L ) ) ,   s h a p e   =   " e l l i p s e " ) ,   
                   f i l l = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) , 
                   a l p h a = 0 . 6 , 
                   q u a n t i t i e s   =   l i s t ( c e x = 2 ) , 
                   l a b e l s   =   l i s t ( c e x = 2 ) , 
                   t i t l e = " " ) 
 
 g < - g r i d . a r r a n g e ( g r o b s   =   l i s t ( p 1 ) ,   t o p   =   t e x t G r o b ( " P e p t i d e s ,   1 %   F D R " ,   g p = g p a r ( f o n t s i z e = 2 5 , f o n t = 8 ) ) ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n - 3 . p d f " , 
               g , 
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 # # # # 
 
 
 # # # #   V e n n   4 
     L F 1 = f f $ p r o t e i n s [ f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
     L F 2 = f f $ p r o t e i n s [ ! f f $ f i l e n a m e % i n % " 2 0 2 5 - 0 3 - 2 1 _ M Y - H R _ L F - D D A _ 2 n g _ n c e 2 8 . m z M L " ] 
 
 T L = u n i q u e ( i n t e r s e c t ( L F 1 , L F 2 ) ) 
 l f 1 < - u n i q u e ( s e t d i f f ( L F 1 , L F 2 ) ) 
 l f 2 < - u n i q u e ( s e t d i f f ( L F 2 , L F 1 ) ) 
 
 p 1 < - p l o t ( e u l e r ( c ( " L F "   =   l e n g t h ( l f 1 ) ,   " T a g "   =   l e n g t h ( l f 2 ) ,   " L F & T a g "   =   l e n g t h ( T L ) ) ,   s h a p e   =   " e l l i p s e " ) ,   
                   f i l l = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " ) , 
                   a l p h a = 0 . 6 , 
                   q u a n t i t i e s   =   l i s t ( c e x = 2 ) , 
                   l a b e l s   =   l i s t ( c e x = 2 ) , 
                   t i t l e = " " ) 
 
 g < - g r i d . a r r a n g e ( g r o b s   =   l i s t ( p 1 ) ,   t o p   =   t e x t G r o b ( " P r o t e i n s ,   1 %   F D R " ,   g p = g p a r ( f o n t s i z e = 2 5 , f o n t = 8 ) ) ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / v e n n - 4 . p d f " , 
               g , 
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 # # # # 