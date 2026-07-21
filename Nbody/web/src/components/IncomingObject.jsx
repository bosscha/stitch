import React, { useEffect, useRef, useState } from 'react';

const IncomingObject = () => {
  const [isVisible, setIsVisible] = useState(false);
  const sectionRef = useRef(null);

  useEffect(() => {
    const observer = new IntersectionObserver(([entry]) => {
      if (entry.isIntersecting) {
        setIsVisible(true);
        observer.disconnect();
      }
    }, { threshold: 0.15 });
    
    if (sectionRef.current) {
      observer.observe(sectionRef.current);
    }
    return () => observer.disconnect();
  }, []);

  return (
    <section ref={sectionRef} id="collisions" className={`glass-section ${isVisible ? 'visible' : ''}`}>
      <h2>Incoming Object Dynamics</h2>
      <p className="section-desc">Testing cluster stability against rogue hyper-massive objects.</p>
      <div className="split-layout">
        <div className="text-content">
          <p>We've implemented specialized simulation forks (<code>rust_nbody_2D_incoming</code> and <code>rust_nbody_3D_incoming</code>) to study the disruption of a perfectly virialized Salpeter-mass cluster.</p>
          <ul className="stats-list">
            <li><strong>Cluster Stars:</strong> 10,000 (King Model, Q=0.5)</li>
            <li><strong>Rogue Mass:</strong> 300 M<sub>&odot;</sub></li>
            <li><strong>Origin Point:</strong> 200 parsecs</li>
            <li><strong>Impact Velocity:</strong> 200 km/s</li>
          </ul>
        </div>
        <div className="visual-placeholder">
          <div className="orbit-ring"></div>
          <div className="cluster"></div>
          <div className="rogue-star"></div>
        </div>
      </div>
    </section>
  );
};

export default IncomingObject;
