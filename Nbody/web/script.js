// Intersection Observer for scroll animations
const observerOptions = {
    root: null,
    rootMargin: '0px',
    threshold: 0.15
};

const observer = new IntersectionObserver((entries, observer) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('visible');
            observer.unobserve(entry.target);
        }
    });
}, observerOptions);

document.querySelectorAll('.glass-section').forEach(section => {
    observer.observe(section);
});

// Interactive Starfield Background
const canvas = document.getElementById('starfield');
const ctx = canvas.getContext('2d');

let width, height;
let particles = [];

function init() {
    width = canvas.width = window.innerWidth;
    height = canvas.height = window.innerHeight;
    particles = [];
    
    const numParticles = Math.floor((width * height) / 8000);
    
    for (let i = 0; i < numParticles; i++) {
        particles.push({
            x: Math.random() * width,
            y: Math.random() * height,
            radius: Math.random() * 1.5,
            color: Math.random() > 0.5 ? '#00F0FF' : (Math.random() > 0.5 ? '#FF0055' : '#FFFFFF'),
            speedX: (Math.random() - 0.5) * 0.5,
            speedY: (Math.random() - 0.5) * 0.5,
            alpha: Math.random()
        });
    }
}

function animate() {
    ctx.clearRect(0, 0, width, height);
    
    particles.forEach(p => {
        p.x += p.speedX;
        p.y += p.speedY;
        
        // Wrap around
        if (p.x < 0) p.x = width;
        if (p.x > width) p.x = 0;
        if (p.y < 0) p.y = height;
        if (p.y > height) p.y = 0;
        
        // Twinkle
        p.alpha += (Math.random() - 0.5) * 0.05;
        if (p.alpha < 0.1) p.alpha = 0.1;
        if (p.alpha > 0.8) p.alpha = 0.8;
        
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2);
        ctx.fillStyle = p.color;
        ctx.globalAlpha = p.alpha;
        ctx.fill();
    });
    
    // Draw connections if close
    ctx.globalAlpha = 0.05;
    ctx.strokeStyle = '#00F0FF';
    ctx.lineWidth = 0.5;
    
    for (let i = 0; i < particles.length; i++) {
        for (let j = i + 1; j < particles.length; j++) {
            const dx = particles[i].x - particles[j].x;
            const dy = particles[i].y - particles[j].y;
            const dist = Math.sqrt(dx*dx + dy*dy);
            
            if (dist < 100) {
                ctx.beginPath();
                ctx.moveTo(particles[i].x, particles[i].y);
                ctx.lineTo(particles[j].x, particles[j].y);
                ctx.stroke();
            }
        }
    }
    
    requestAnimationFrame(animate);
}

window.addEventListener('resize', init);
init();
animate();
