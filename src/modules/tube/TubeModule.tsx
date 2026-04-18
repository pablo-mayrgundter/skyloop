import { Canvas } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import LoopTube from '../../three/loop/LoopTube'

const TubeModule = () => {
  return (
    <div className="canvas-wrapper" style={{ height: 320 }}>
      <Canvas camera={{ position: [0.2, 0.3, 0.8] }}>
        <ambientLight intensity={0.6} />
        <LoopTube />
        <OrbitControls />
      </Canvas>
      <p>Examining tube geometry, evacuation strategy, and structural thickness.</p>
    </div>
  )
}

export default TubeModule
