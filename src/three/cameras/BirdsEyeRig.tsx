import {OrbitControls, PerspectiveCamera} from '@react-three/drei'
import {useRef} from 'react'
import {PerspectiveCamera as PerspectiveCameraImpl} from 'three'

const BirdsEyeRig = () => {
  const camRef = useRef<PerspectiveCameraImpl>(null)
  return (
    <>
      <PerspectiveCamera ref={camRef} makeDefault position={[0, 0.15, 0.2]} near={0.001} far={200} fov={50} />
      <OrbitControls enablePan enableZoom maxDistance={1} />
    </>
  )
}

export default BirdsEyeRig
