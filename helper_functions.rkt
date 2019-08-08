#lang racket
(require "declarations.rkt")

(provide buildTree calcForces moveparticles)

(define (cal-lx area)
  (bbox-llx area))

(define (cal-ux area)
  (bbox-rux area))

(define (cal-ly area)
  (bbox-lly area))

(define (cal-uy area)
  (bbox-ruy area))


(define (system-mass particles)
  (foldr + 0 (lc (particle-mass p) : p <- particles))) 

(define (centre-of-mass particles)     
  (vec (/ (foldr (lambda(x y) (+ (* (particle-mass x) (vec-x (particle-posn x))) y)) 0 particles) (system-mass particles))
       (/ (foldr (lambda(x y) (+ (* (particle-mass x) (vec-y (particle-posn x))) y)) 0 particles) (system-mass particles))))

(define (buildTree initialArea particles)
 (let*  ((lx ( cal-lx initialArea ))
         (ly ( cal-ly initialArea))
         (ux (cal-ux initialArea))
         (uy (cal-uy initialArea))
         (px (/ (+ lx ux) 2))
         (py (/ (+ ly uy) 2))
         (p1area (bbox lx py px uy))
         (p2area (bbox px py ux uy))
         (p3area (bbox lx ly px py))
         (p4area (bbox px ly ux py))
         (p1-particles (lc p : p <- particles @ (and (< (vec-x (particle-posn p)) px) (>= (vec-y (particle-posn p)) py))))    
         (p2-particles (lc p : p <- particles @ (and (>= (vec-x (particle-posn p)) px) (>= (vec-y (particle-posn p)) py))))        
         (p3-particles (lc p : p <- particles @ (and (< (vec-x (particle-posn p)) px) (< (vec-y (particle-posn p)) py))))
         (p4-particles (lc p : p <- particles @ (and (>= (vec-x (particle-posn p)) px) (< (vec-y (particle-posn p)) py)))))
   (cond [(null? particles) '()]
         [(singleton particles)  (gnode (particle-mass (car particles)) (particle-posn (car particles)) '())]
         [#t (gnode (system-mass particles) (centre-of-mass particles)
          (remove* '(())  (list (buildTree p1area p1-particles) (buildTree p2area p2-particles) (buildTree p3area p3-particles) (buildTree p4area p4-particles))))])))

(define (is-near? single-particle tree s)
      (if (<= (distance (particle-posn single-particle) (gnode-posn tree)) (* theta s)) #t
          #f))

(define (distance a b)
  (sqrt (+ (expt (- (vec-x a) (vec-x b)) 2) (expt (- (vec-y a) (vec-y b)) 2))))

(define (add-vectors a b)
  (vec (+ (vec-x a) (vec-x b)) (+ (vec-y a) (vec-y b))))

(define (mult-vector a t)
  (vec (* (vec-x a) t) (* (vec-y a) t)))

(define (single-particle-force particle1 particle2-mass particle2-posn)
  (let* [(force-x-comp (/ (* g (particle-mass particle1) particle2-mass (- (vec-x particle2-posn) (vec-x (particle-posn particle1))))
                   (expt (distance (particle-posn particle1) particle2-posn) 3)))
         (force-y-comp (/ (* g (particle-mass particle1) particle2-mass (- (vec-y particle2-posn) (vec-y (particle-posn particle1))))
                   (expt (distance (particle-posn particle1) particle2-posn) 3)))]
    (vec force-x-comp force-y-comp)))

(define (calcForces initialArea tree particles)
  (let ((side-of-square (abs (- (bbox-rux initialArea) (bbox-llx initialArea)))))
    (define (total-force-sp p tree side-length)
      (cond [(null? (gnode-subtrees tree)) (if (and (equal? (particle-posn p) (gnode-posn tree)) (=  (particle-mass p) (gnode-mass tree))) (vec 0 0)
                                             (single-particle-force p (gnode-mass tree) (gnode-posn tree)))]
            [#t (if (is-near? p tree side-length) (foldr (lambda(x y) (add-vectors (total-force-sp p x (/ side-length 2)) y)) (vec 0 0) (gnode-subtrees tree))
                    (single-particle-force p (gnode-mass tree) (gnode-posn tree)))]))
    (map (lambda(x) (total-force-sp x tree side-of-square)) particles)))

(define (moveparticles particles forces)
  (if (null? particles) '()
  (let* ((p (car particles))
         (f (car forces)))
  (define (move-single-particle  particle1 force)
    (let* ((acceleration (vec (/ (vec-x force) (particle-mass particle1)) (/ (vec-y force) (particle-mass particle1))))
           (displacement (add-vectors (mult-vector (particle-velocity particle1) timeslice) (mult-vector acceleration (* 0.5 (expt timeslice 2)))))
           (final-posn (add-vectors displacement (particle-posn particle1)))
           (final-velocity (add-vectors (particle-velocity particle1) (mult-vector acceleration timeslice))))
      (particle (particle-mass particle1) final-posn final-velocity)))
    (cons (move-single-particle p f) (moveparticles (cdr particles) (cdr forces))))))



