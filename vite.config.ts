import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react-swc'

// Determine base path: use repo name for GitHub Pages, '/' for local dev
const getBasePath = () => {
  if (process.env.GITHUB_REPOSITORY) {
    const repoName = process.env.GITHUB_REPOSITORY.split('/')[1]
    return `/${repoName}/`
  }
  return '/'
}

export default defineConfig({
  base: getBasePath(),
  plugins: [
    react(),
    {
      name: 'fix-three-js-compatibility',
      enforce: 'pre',
      transform(code, id) {
        // Skip node_modules except the ones we need to fix
        if (!id.includes('node_modules')) {
          return null
        }
        
        // Fix deprecated BufferGeometry imports in troika packages and drei
        if (id.includes('troika-three-utils') || id.includes('troika-three-text') || id.includes('@react-three/drei') || id.includes('three-stdlib')) {
          let transformed = code
            .replace(/CylinderBufferGeometry/g, 'CylinderGeometry')
            .replace(/BoxBufferGeometry/g, 'BoxGeometry')
            .replace(/PlaneBufferGeometry/g, 'PlaneGeometry')
            .replace(/SphereBufferGeometry/g, 'SphereGeometry')
            .replace(/RingBufferGeometry/g, 'RingGeometry')
            .replace(/TorusBufferGeometry/g, 'TorusGeometry')
          
          // Fix deprecated encoding constants (but skip three-stdlib as it has its own declarations)
          if (id.includes('@react-three/drei') && !id.includes('three-stdlib')) {
            // Replace in import statements from 'three'
            transformed = transformed.replace(/import\s*\{([^}]*)\}\s*from\s*['"]three['"]/g, (match, imports) => {
              const newImports = imports
                .replace(/sRGBEncoding/g, 'SRGBColorSpace')
                .replace(/LinearEncoding/g, 'LinearSRGBColorSpace')
                .replace(/LuminanceFormat/g, 'RedFormat')
              return match.replace(imports, newImports)
            })
            // Replace usage (but not declarations)
            transformed = transformed.replace(/\bsRGBEncoding\b/g, 'SRGBColorSpace')
            transformed = transformed.replace(/\bLinearEncoding\b/g, 'LinearSRGBColorSpace')
            transformed = transformed.replace(/\bLuminanceFormat\b/g, 'RedFormat')
          }
          
          // Fix three-stdlib separately - remove LuminanceFormat from imports and replace usage with RedFormat
          if (id.includes('three-stdlib')) {
            // Remove LuminanceFormat from import statements from 'three' (RedFormat is already imported)
            transformed = transformed.replace(/import\s*\{([^}]*)\}\s*from\s*['"]three['"]/g, (match, imports) => {
              if (imports.includes('LuminanceFormat')) {
                const newImports = imports
                  .split(',')
                  .map(i => i.trim())
                  .filter(i => i !== 'LuminanceFormat')
                  .join(', ')
                return match.replace(imports, newImports)
              }
              return match
            })
            // Replace usage of LuminanceFormat with RedFormat
            transformed = transformed.replace(/\bLuminanceFormat\b/g, 'RedFormat')
          }
          
          return {
            code: transformed,
            map: null
          }
        }
      }
    }
  ],
  server: {
    port: 5173
  },
  optimizeDeps: {
    esbuildOptions: {
      plugins: [
        {
          name: 'fix-three-js-compatibility-esbuild',
          setup(build) {
            // Use a more specific filter for the files we need to transform
            const filter = /(troika-three-utils|troika-three-text|@react-three\/drei|three-stdlib).*\.(js|mjs|esm\.js)$/
            build.onLoad({ filter }, async (args) => {
              const fs = await import('fs')
              let contents = await fs.promises.readFile(args.path, 'utf8')
              
              // Fix deprecated BufferGeometry imports
              contents = contents
                .replace(/CylinderBufferGeometry/g, 'CylinderGeometry')
                .replace(/BoxBufferGeometry/g, 'BoxGeometry')
                .replace(/PlaneBufferGeometry/g, 'PlaneGeometry')
                .replace(/SphereBufferGeometry/g, 'SphereGeometry')
                .replace(/RingBufferGeometry/g, 'RingGeometry')
                .replace(/TorusBufferGeometry/g, 'TorusGeometry')
              
              // Fix deprecated encoding constants in drei
              if (args.path.includes('@react-three/drei') && !args.path.includes('three-stdlib')) {
                contents = contents.replace(/import\s*\{([^}]*)\}\s*from\s*['"]three['"]/g, (match, imports) => {
                  const newImports = imports
                    .replace(/sRGBEncoding/g, 'SRGBColorSpace')
                    .replace(/LinearEncoding/g, 'LinearSRGBColorSpace')
                    .replace(/LuminanceFormat/g, 'RedFormat')
                  return match.replace(imports, newImports)
                })
                contents = contents.replace(/\bsRGBEncoding\b/g, 'SRGBColorSpace')
                contents = contents.replace(/\bLinearEncoding\b/g, 'LinearSRGBColorSpace')
                contents = contents.replace(/\bLuminanceFormat\b/g, 'RedFormat')
              }
              
              // Fix three-stdlib separately
              if (args.path.includes('three-stdlib')) {
                contents = contents.replace(/import\s*\{([^}]*)\}\s*from\s*['"]three['"]/g, (match, imports) => {
                  if (imports.includes('LuminanceFormat')) {
                    const newImports = imports
                      .split(',')
                      .map(i => i.trim())
                      .filter(i => i !== 'LuminanceFormat')
                      .join(', ')
                    return match.replace(imports, newImports)
                  }
                  return match
                })
                contents = contents.replace(/\bLuminanceFormat\b/g, 'RedFormat')
              }
              
              return { contents, loader: 'js' }
            })
          }
        }
      ]
    }
  },
  test: {
    globals: true,
    environment: 'jsdom',
    exclude: ['**/node_modules/**', '**/dist/**', '**/e2e/**', '**/*.e2e.{test,spec}.{js,mjs,cjs,ts,mts,cts,jsx,tsx}']
  }
})
