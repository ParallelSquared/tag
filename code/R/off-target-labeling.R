#   E v a l u a t i n g   o f f - t a r g e t   l a b e l i n g   w i t h   P S M t a g 
 
 l i b r a r y ( r e s h a p e 2 ) 
 
 e v < - r e a d . d e l i m ( " / V o l u m e s / L a b / S u b m i s s i o n s / t a g _ p r e p r i n t _ v 1 / O f f - t a r g e t - l a b e l i n g / c o m b i n e d / t x t / e v i d e n c e . t x t " ) 
 
 d f < - m e l t ( e v [ , 1 7 : 2 0 ] ) 
 
 d f $ v a r i a b l e < - g s u b ( " T 6 _ " , " " , d f $ v a r i a b l e ) 
 
 p d f < -   d f   % > %   g r o u p _ b y ( v a r i a b l e )   % > %   s u m m a r i s e ( p c t = m e a n ( v a l u e ) ) 
 
 g g p l o t ( p d f , a e s ( x = v a r i a b l e , y = 1 0 0 * p c t ) )   +   
     g e o m _ b a r ( s t a t = " i d e n t i t y " )   +   
     t h e m e _ t a g ( )   +   
     x l a b ( " \ n P o s s i b l e   l a b e l i n g   s i t e " )   +   
     y l a b ( " %   l a b e l e d \ n " )   +   
     g g t i t l e ( " H u m a n   p r o t e o m e ,   4   n g " ,   s u b t i t l e = " 1 3 , 8 6 5   u n i q u e   p r e c u r s o r s \ n 3   t e c h n i c a l   r e p l i c a t e s " )   
 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / o f f - t a r g e t - H S T Y . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 
 e v < - r e a d . d e l i m ( " / V o l u m e s / L a b / H S / 2 0 2 5 0 4 3 0 _ t 6 - o p e n / p t m - s h e p h e r d - o u t p u t / g l o b a l . m o d s u m m a r y . t s v " ) 
 
 e v $ m s _ r o u n d e d < - r o u n d ( e v $ M a s s . S h i f t , 1 ) 
 
 p d f < - e v   % > %   g r o u p _ b y ( m s _ r o u n d e d )   % > %   s u m m a r i s e ( f r e q   =   s u m ( d a t a s e t 0 1 _ p e r c e n t _ P S M s ) ) 
 
 g g p l o t ( p d f , a e s ( x = a s . f a c t o r ( m s _ r o u n d e d ) , y = f r e q ) )   +   
     g e o m _ b a r ( s t a t = " i d e n t i t y " )   +   
     t h e m e _ t a g ( )   +   
     x l a b ( " \ n M o d i f i c a t i o n   m a s s " )   +   
     y l a b ( " %   o f   o b s e r v e d   P S M s \ n " )   +   
     g g t i t l e ( " O p e n   s e a r c h ,   - 5 0   t o   + 1 3 0 0 D a " ,   s u b t i t l e = " H u m a n   p r o t e o m e ,   4   n g \ n 3   t e c h n i c a l   r e p l i c a t e s " )   
 
 
 g g s a v e ( " / U s e r s / h s / L i b r a r y / C l o u d S t o r a g e / G o o g l e D r i v e - h s p e c h t @ p a r a l l e l s q . o r g / S h a r e d   d r i v e s / P T I   S h a r e d   D r i v e / P r o j e c t s / T e a m T a g / t a g 6 _ p u b l i c a t i o n / f i g u r e s / o p e n - s e a r c h . p d f " ,   
               d e v i c e = " p d f " ,   
               w i d t h   =   i m a g e _ i n _ w ,   
               h e i g h t   =   i m a g e _ i n _ h , 
               u n i t s = " i n " ) 
 
 
 