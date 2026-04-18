import { Canvas } from '@react-three/fiber'
import LoopSections from '../../three/loop/LoopSections'
import Globe from '../../three/Globe'

const SectionsModule = () => (
  <div className="panel">
    <h3>Sections Lab</h3>
    <p>Study buoyant vacuum sections along the corridor.</p>
    <div className="canvas-wrapper" style={{ height: '320px' }}>
      <Canvas camera={{ position: [0.05, 0.05, 0.15] }}>
        <ambientLight intensity={0.5} />
        <Globe />
        <LoopSections count={10} elevationKm={15} />
      </Canvas>
    </div>
  </div>
)

export default SectionsModule
