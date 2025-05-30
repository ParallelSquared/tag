#   E v a l u a t i n g   a p p l i c a t i o n   o f   T a g   t o   h u m a n   p r o t e o m e   -   S C P   d a t a 
 
 s o u r c e ( " p a r a m s . R " ) 
 
 #   R e a d   i n   l i s t   o f   y e a s t   p r o t e i n s   
 y e a s t < - r e a d . d e l i m ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / f a s t a / u n i p r o t k b _ r e v i e w e d _ t r u e _ A N D _ m o d e l _ o r g a n _ 2 0 2 5 _ 0 4 _ 0 4 . t s v " ) 
 
 #   R e a d   i n   M S   d a t a 
 l f < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - l a b e l - f r e e / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - s i n g l e - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 5 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - 5 - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
   
 #   R e m o v e   y e a s t   p r o t e i n s   f r o m   d a t a 
 l f < - l f [ ! s t r _ s p l i t ( l f $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 < - t 6 [ ! s t r _ s p l i t ( t 6 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 5 < - t 6 5 [ ! s t r _ s p l i t ( t 6 5 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 
 #   D u p l i c a t e   d a t a 
 l f 2 0 0 < - l f [ g r e p l ( " 2 0 0 p g " , l f $ R u n ) , ] 
 t 2 0 0 < - t 6 # [ g r e p l ( " 2 0 0 p g _ 2 0 w i n " , t 6 $ R u n ) , ] 
 t 5 _ 2 0 0 < - t 6 5 # t 6 [ g r e p l ( " 2 0 0 p g _ 5 p l e x _ 2 0 w i n " , t 6 $ R u n ) , ] 
 
 #   G r a b   p r e c u r s o r s   c o n t a i n i n g   R 
 l f 2 0 0 < - l f 2 0 0 [ g r e p l ( " R " , l f 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 t 2 0 0 < - t 2 0 0 [ g r e p l ( " R " , t 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 t 5 _ 2 0 0 < - t 5 _ 2 0 0 [ g r e p l ( " R " , t 5 _ 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 
 #   D e f i n e   p r e c u r s o r   m o d i f i e d   s e q u e n c e   +   c h a r g e   i n   o n e   s t r i n g 
 l f 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( l f 2 0 0 $ M o d i f i e d . S e q u e n c e ,   l f 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 5 _ 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 
 #   R e m o v e   d u p l i c a t e   e n t r i e s   p e r   r e p l i c a t e ,   i f   a n y 
 l f 2 0 0 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 5 _ 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 
 #   C o u n t   p r e c u r s o r s 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 0 0 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 2 0 0 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 5 _ 2 0 0 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 3 , 2 ) , 1 , 2 , 3 ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , " T a g - 5 p l e x " , " T a g - 5 p l e x " , " T a g - 5 p l e x " ) 
 
 
 d f x < - d f 
 
 #   R e m o v e   d u p l i c a t e   p r o t e i n   e n t r i e s   p e r   r e p l i c a t e ,   i f   a n y 
 l f 2 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 3 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 2 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 
 #   C o u n t   p r o t e i n s   - -   n o   F D R   f i l t e r i n g ,   t h i s   i s   j u s t   f o r   i l l u s t r a t i n g   e f f e c t   o f 
 #   o n e   v .   t w o   t a g s ,   n o t   t o t a l   c o u n t s   a c r o s s   s a m p l e   t y p e s 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 6 2 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 6 3 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 2 , 3 ) , 1 , 2 , 3 ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , " T a g - 5 p l e x " , " T a g - 5 p l e x " , " T a g - 5 p l e x " ) 
 
 d f 3 < - r b i n d ( d f x , d f ) 
 d f 3 $ T y p e < - c ( r e p ( " P r e c u r s o r s " , 9 ) , r e p ( " P r o t e i n   I D s " , 9 ) ) 
 
 d f 4 < - d f 3   % > %   g r o u p _ b y ( L a b e l ,   T y p e )   % > %   s u m m a r i s e ( s d   =   s d ( F r e q ) ,   F r e q m   =   m e a n ( F r e q ) )   % > %   u n g r o u p ( ) 
 
 # # # #   R e p e a t   a b o v e   f o r   K   c o n t a i n i n g   p r e c u r s o r s   a n d   p r o t e i n s 
 y e a s t < - r e a d . d e l i m ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / f a s t a / u n i p r o t k b _ r e v i e w e d _ t r u e _ A N D _ m o d e l _ o r g a n _ 2 0 2 5 _ 0 4 _ 0 4 . t s v " ) 
 l f < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - l a b e l - f r e e / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - s i n g l e - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 5 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - 5 - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 
 l f < - l f [ ! s t r _ s p l i t ( l f $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 < - t 6 [ ! s t r _ s p l i t ( t 6 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 5 < - t 6 5 [ ! s t r _ s p l i t ( t 6 5 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 
 l f 2 0 0 < - l f [ g r e p l ( " 2 0 0 p g " , l f $ R u n ) , ] 
 t 2 0 0 < - t 6 # [ g r e p l ( " 2 0 0 p g _ 2 0 w i n " , t 6 $ R u n ) , ] 
 t 5 _ 2 0 0 < - t 6 5 # t 6 [ g r e p l ( " 2 0 0 p g _ 5 p l e x _ 2 0 w i n " , t 6 $ R u n ) , ] 
 
 l f 2 0 0 < - l f 2 0 0 [ g r e p l ( " K " , l f 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 t 2 0 0 < - t 2 0 0 [ g r e p l ( " K " , t 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 t 5 _ 2 0 0 < - t 5 _ 2 0 0 [ g r e p l ( " K " , t 5 _ 2 0 0 $ S t r i p p e d . S e q u e n c e ) , ] 
 
 l f 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( l f 2 0 0 $ M o d i f i e d . S e q u e n c e ,   l f 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 5 _ 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 
 l f 2 0 0 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 5 _ 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 0 0 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 2 0 0 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 5 _ 2 0 0 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 3 , 2 ) , 1 , 2 , 3 ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , " T a g - 5 p l e x " , " T a g - 5 p l e x " , " T a g - 5 p l e x " ) 
 
 d f x < - d f 
 
 l f 2 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 3 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 2 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 6 2 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 6 3 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 3 , 2 ) , 1 , 2 , 3 ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , " T a g - 5 p l e x " , " T a g - 5 p l e x " , " T a g - 5 p l e x " ) 
 
 d f 3 < - r b i n d ( d f x , d f ) 
 d f 3 $ T y p e < - c ( r e p ( " P r e c u r s o r s " , 9 ) , r e p ( " P r o t e i n   I D s " , 9 ) ) 
 d f 4 k < - d f 3   % > %   g r o u p _ b y ( L a b e l ,   T y p e )   % > %   s u m m a r i s e ( s d   =   s d ( F r e q ) ,   F r e q m   =   m e a n ( F r e q ) )   % > %   u n g r o u p ( ) 
 d f 4 $ ` C - t e r m   a m i n o   a c i d `   < - " R " 
 d f 4 k $ ` C - t e r m   a m i n o   a c i d `   < - " K " 
 d f 4 r k < - r b i n d ( d f 4 , d f 4 k ) ;   d f 4 r k 
 d f 4 r k $ ` C - t e r m   a m i n o   a c i d ` < - f a c t o r ( d f 4 r k $ ` C - t e r m   a m i n o   a c i d ` ,   l e v e l s   =   c ( " R " , " K " ) ) 
 
 
 # # #   V i s u a l i z e !   
 
 g g p l o t ( d f 4 r k , a e s ( x = L a b e l , y = F r e q m , f i l l = L a b e l ,   a l p h a = ` C - t e r m   a m i n o   a c i d ` ) )   + 
     g e o m _ b a r ( s t a t = " i d e n t i t y " ,   p o s i t i o n = " d o d g e " )   + 
     s c a l e _ a l p h a _ d i s c r e t e ( r a n g e = c ( 0 . 7 , 0 . 3 ) )   + 
     g e o m _ e r r o r b a r ( a e s ( y m i n = F r e q m - s d ,   y m a x = F r e q m + s d ) ,   w i d t h = . 2 , 
                                 p o s i t i o n = p o s i t i o n _ d o d g e ( . 9 ) )   + 
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " , " # f f 7 f 1 e " ) )   +   
     t h e m e _ t a g ( )   +   
     t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     # r r e m o v e ( " l e g e n d " )   +   
     g u i d e s ( f i l l = " n o n e " )   + 
     x l a b ( " " )   + 
     y l a b ( " C o u n t ,   1   % F D R \ n " )   +   
     f a c e t _ w r a p ( " T y p e " ,   s c a l e s   =   " f r e e " )   +   
     g g t i t l e ( " t i m s T O F   S C P ,   D I A - N N " ) #   +   
     #   t h e m e ( 
     #       p l o t . t i t l e . p o s i t i o n   =   " p l o t " ,   
     #       p l o t . c a p t i o n . p o s i t i o n   =     " p l o t " 
     #   ) 
 
 #   S a v e   t o   G o o g l e   D r i v e 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / d i a _ i d _ s c p . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h * 1 . 5 , 
               u n i t s = " i n " ) 
 
 # # # # # # # # # # # # #   T o t a l   I D s ,   s t a c k e d   
 
 y e a s t < - r e a d . d e l i m ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / f a s t a / u n i p r o t k b _ r e v i e w e d _ t r u e _ A N D _ m o d e l _ o r g a n _ 2 0 2 5 _ 0 4 _ 0 4 . t s v " ) 
 l f < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - l a b e l - f r e e / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - s i n g l e - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 t 6 5 < - r e a d _ p a r q u e t ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 5 2 0 _ P S M t a g _ s h a r e / s i n g l e c e l l d a t a - s c p / d i a n n - r e s u l t s - 5 - p l e x / r e p o r t - f i r s t - p a s s . p a r q u e t " ) 
 
 l f < - l f [ ! s t r _ s p l i t ( l f $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 < - t 6 [ ! s t r _ s p l i t ( t 6 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 t 6 5 < - t 6 5 [ ! s t r _ s p l i t ( t 6 5 $ P r o t e i n . G r o u p , " ; " )   % i n %   y e a s t $ E n t r y ,   ] 
 
 l f 2 0 0 < - l f [ g r e p l ( " 2 0 0 p g " , l f $ R u n ) , ] 
 t 2 0 0 < - t 6 # [ g r e p l ( " 2 0 0 p g _ 2 0 w i n " , t 6 $ R u n ) , ] 
 t 5 _ 2 0 0 < - t 6 5 # t 6 [ g r e p l ( " 2 0 0 p g _ 5 p l e x _ 2 0 w i n " , t 6 $ R u n ) , ] 
 
 t 5 _ 2 0 0 $ R u n < - p a s t e 0 ( t 5 _ 2 0 0 $ R u n ,   t 5 _ 2 0 0 $ C h a n n e l ) 
 
 l f 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( l f 2 0 0 $ M o d i f i e d . S e q u e n c e ,   l f 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e < - p a s t e 0 ( t 5 _ 2 0 0 $ M o d i f i e d . S e q u e n c e ,   t 5 _ 2 0 0 $ P r e c u r s o r . C h a r g e ) 
 
 l f 2 0 0 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 t 5 _ 2 0 0 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " M o d i f i e d . S e q u e n c e " ) ) 
 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 0 0 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 2 0 0 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 5 _ 2 0 0 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 3 , 2 ) , r e p ( 1 : 5 , 3 ) ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , r e p ( p a s t e 0 ( " T a g - 5 p l e x " , " _ " , 1 : 5 ) , 3 ) ) 
 
 d f x < - d f 
 
 #   F i l t e r   t o   1 %   g l o b a l   F D R   ( p e r   e x p e r i m e n t ) 
 l f 2 0 0 < - l f 2 0 0   % > %   g r o u p _ b y ( P r o t e i n . G r o u p )   % > %   d p l y r : : m u t a t e ( " m i n _ P G Q v a l "   =   m i n ( P G . Q . V a l u e , n a . r m = T ) )   % > %   u n g r o u p ( ) 
 t 5 _ 2 0 0 < - t 5 _ 2 0 0   % > %   g r o u p _ b y ( P r o t e i n . G r o u p )   % > %   d p l y r : : m u t a t e ( " m i n _ P G Q v a l "   =   m i n ( P G . Q . V a l u e , n a . r m = T ) )   % > %   u n g r o u p ( ) 
 t 2 0 0 < - t 2 0 0   % > %   g r o u p _ b y ( P r o t e i n . G r o u p )   % > %   d p l y r : : m u t a t e ( " m i n _ P G Q v a l "   =   m i n ( P G . Q . V a l u e , n a . r m = T ) )   % > %   u n g r o u p ( ) 
 
 l f 2 0 0 < - l f 2 0 0 [ l f 2 0 0 $ m i n _ P G Q v a l < 0 . 0 1 , ] 
 t 5 _ 2 0 0 < - t 5 _ 2 0 0 [ t 5 _ 2 0 0 $ m i n _ P G Q v a l < 0 . 0 1 , ] 
 t 2 0 0 < - t 2 0 0 [ t 2 0 0 $ m i n _ P G Q v a l < 0 . 0 1 , ] 
 
 l f 2 < - r e m o v e . d u p l i c a t e s ( l f 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 3 < - r e m o v e . d u p l i c a t e s ( t 5 _ 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 t 6 2 < - r e m o v e . d u p l i c a t e s ( t 2 0 0 , c ( " R u n " , " P r o t e i n . I d s " ) ) 
 
 d f 1 < - d a t a . f r a m e ( t a b l e ( l f 2 $ R u n ) ) 
 d f 2 < - d a t a . f r a m e ( t a b l e ( t 6 2 $ R u n ) ) 
 d f 3 < - d a t a . f r a m e ( t a b l e ( t 6 3 $ R u n ) ) 
 
 d f < - r b i n d ( d f 1 , d f 2 ) 
 d f < - r b i n d ( d f , d f 3 ) 
 d f $ R e p < - a s . c h a r a c t e r ( c ( r e p ( 1 : 3 , 2 ) , r e p ( 1 : 5 , 3 ) ) ) 
 d f $ L a b e l < - c ( " L a b e l - f r e e " , " L a b e l - f r e e " , " L a b e l - f r e e " , " T a g " , " T a g " , " T a g " , r e p ( p a s t e 0 ( " T a g - 5 p l e x " , " _ " , 1 : 5 ) , 3 ) ) 
 
 d f 3 < - r b i n d ( d f x , d f ) 
 d f 3 $ T y p e < - c ( r e p ( " P r e c u r s o r s " , 2 1 ) , r e p ( " P r o t e i n   I D s " , 2 1 ) ) 
 
 d f 4 < - d f 3   % > %   g r o u p _ b y ( L a b e l ,   T y p e )   % > %   s u m m a r i s e ( s d   =   s d ( F r e q ) ,   F r e q m   =   m e a n ( F r e q ) )   % > %   u n g r o u p ( ) 
 
 d f 4 $ L T < - c ( r e p ( " L a b e l - f r e e " , 2 ) , r e p ( " T a g " , 2 ) , r e p ( " T a g - 5 p l e x " , 1 0 ) ) 
 d f 4 $ a l p < - c ( 1 , 1 , 1 , 1 , 2 , 2 , 3 , 3 , 4 , 4 , 5 , 5 , 6 , 6 ) 
 
 e r r o r _ b a r s   < -   d f 4   % > % 
     a r r a n g e ( L T ,   d e s c ( L a b e l ) )   % > % 
     g r o u p _ b y ( L T , T y p e )   % > % 
     m u t a t e ( m e a n _ n e w   =   c u m s u m ( F r e q m ) )   % > % 
     u n g r o u p ( ) 
 
 c u s t o m _ y   < -   l i s t ( 
     s c a l e _ y _ c o n t i n u o u s ( l i m i t s   =   c ( 0 ,   7 0 0 0 0 ) , l a b e l s   =   s c a l e s : : c o m m a , b r e a k s   =   s c a l e s : : p r e t t y _ b r e a k s ( n   =   6 ) ) , 
     s c a l e _ y _ c o n t i n u o u s ( l i m i t s   =   c ( 0 ,   1 0 0 0 0 ) , l a b e l s   =   s c a l e s : : c o m m a , b r e a k s   =   s c a l e s : : p r e t t y _ b r e a k s ( n   =   5 ) ) 
 ) 
 
 
 g g p l o t ( d f 4 , a e s ( x = L T , y = F r e q m , f i l l = L a b e l , a l p h a = a l p ) )   + 
     g e o m _ b a r ( s t a t = " i d e n t i t y " , a e s ( f i l l = L a b e l , a l p h a = a l p ) )   + 
     s c a l e _ a l p h a _ c o n t i n u o u s ( r a n g e = c ( 1 , 0 . 4 ) )   + 
     g e o m _ e r r o r b a r ( d a t a = e r r o r _ b a r s ,   a e s ( x = L T ,   y m i n = m e a n _ n e w - s d ,   y m a x = m e a n _ n e w + s d ) ,   w i d t h = . 2 )   +   
     s c a l e _ f i l l _ m a n u a l ( v a l u e s = c ( " # 1 f 7 7 b 4 " , " # f f 7 f 0 e " , r e p ( " # f f 7 f 0 e " , 5 ) ) )   +   
     s c a l e _ c o l o r _ m a n u a l ( v a l u e s = r e p ( " g r e y 3 0 " , 7 ) )   +   
     t h e m e _ t a g ( )   +   
     t h e m e ( a x i s . t e x t . x   =   e l e m e n t _ t e x t ( a n g l e   =   4 5 ,   v j u s t   =   1 ,   h j u s t = 1 ) )   +   
     r r e m o v e ( " l e g e n d " )   +   
     x l a b ( " " )   + 
     y l a b ( " C o u n t ,   1   % F D R \ n " )   +   
     # s c a l e _ f i l l _ g r e y ( )   +   
     f a c e t _ w r a p ( " T y p e " ,   s c a l e s   =   " f r e e _ y " ) +   
     f a c e t t e d _ p o s _ s c a l e s ( y   =   c u s t o m _ y )   +   
     g g t i t l e ( " 1 5   m i n u t e   L C   g r a d i e n t " ,   s u b t i t l e   =   " t i m s T O F   S C P ,   D I A - N N " )   + 
     s c a l e _ y _ c o n t i n u o u s ( b r e a k s   =   s c a l e s : : p r e t t y _ b r e a k s ( n   =   7 ) )   
 #   t h e m e ( 
 #       p l o t . t i t l e . p o s i t i o n   =   " p l o t " ,   
 #       p l o t . c a p t i o n . p o s i t i o n   =     " p l o t " 
 #   ) 
 
 
 d f 4   % > %   g r o u p _ b y ( L T , T y p e )   % > %   s u m m a r i s e ( s u m ( F r e q m ) ) 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / s t a c k e d _ i d s _ s c p . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h * 1 . 5 , 
               u n i t s = " i n " ) 
 