import { Canvas } from '@react-three/fiber'
import LoopGliders from '../../three/loop/LoopGliders'
import Globe from '../../three/Globe'

const GlidersModule = () => (
  <div className="panel">
    <h3>Gliders Lab</h3>
    <p>Stubbed glider dynamics driven by rapier hooks.</p>
    <div className="canvas-wrapper" style={{ height: '320px' }}>
      <Canvas camera={{ position: [0.05, 0.06, 0.16] }}>
        <ambientLight intensity={0.8} />
        <Globe />
        <LoopGliders count={4} />
      </Canvas>
    </div>
  </div>
)

export default GlidersModule
